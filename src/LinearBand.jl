module LinearBand

export linearband

using StaticArrays

# Structure to store segment information.
# A Segment has a first and last index into the t and y arrays indicating the first and last
# points in this segment.  It also stores the segment's slope and height, plus an array
# of three indices for the pivot points.  These indices are relative to the segment, so if the 
# first pivot is 1, this refers to the point t[segment.first].  Finally there is an array of 
# three integers that specify the type for each of the pivots.  0 = undefined, 1 = upper, and
# -1 = lower (meaning the pivot is on the lower boundary of the band).
struct Segment
	first::Int64  # index into t,y arrays of first point in this segment
	last::Int64   # index into t,y arrays of last point in this segment
	splittable::Bool
	slope::Float64 
	height::Float64
	pivot::MVector{3,Int64}  # array of size 3, indices relative to this segment
	pivottype::MVector{3,Int64}  # array of size 3, 0=unknown, 1=upper, -1=lower
end

# Structure to store solution information.
# A Soln has segsort, and tjoin arrays, and the height and tmin values.
# segsort keeps track of all the segments sorted by time, segment[segsort[i]] is the i'th 
# segment from the left for this solution.
# tjoin is an array (with length(segsort)+1 entries, that stores the join points between 
# adjacent unswallowed segments.
# tjoin[i] is the intersection of the midlines of the i'th segment with the first
# unswallowed segment to the left.
# yjoin is an array that stores the y-values of the midline of the segments at the points tjoin
# By definition the join point for the first segment is t[1], and tjoin[nseg+1] is t[n].
struct Soln
	segsort::Vector{Int64}
	tjoin::Vector{Float64}
	yjoin::Vector{Float64}
	height::Float64
	tmin::Float64
end

include("sandwich.jl")
include("split.jl")
include("join.jl")
using Plots

@enum Mode reduction recovery smoothing

"""
	(band,height,tmin) = linearband(t, y, tmin_constraint; minheightfrac=0.8, variableheight=false)

Compute a piecewise linear band around the data (`t`,`y`).
The band is a near minimax continuous piecewise linear band enclosing the data.  No
segment of this band whose slope does not lie between its neighbour's slopes will have
midline length less than `tmin_constraint`.  The returned value `tmin` is the minimum such
length in the actual solution and the returned value `height` is the maximal height of the
band.  The output `band` is a tuple (LB,UB).  The band's lower bound is defined by
connecting the points in the mx2 matrix LB (first column is the t values, second column
is the y values) and the upper bound by the nx2 matrix UB.   The band may have constant
height, or each linear segment may have its own height (if `varableheight=true`).  In the
case that a constant height band is requested, the t values in both LB and UB will be
identical; in the variable height case it is possible that m is different from n.  In
either case, the first and last values of the first column of LB and UB will be the same
as the first and last values of the input `t`.  

The variable `minheightfrac` controls how much splitting of segments is done after the
largest height segment that cannot be split without causing a violation of the constraints
is found.  Segments whose height is above `minheightfrac*height` will be split.  In the
case that `variableheight=false`, all segments are expanded vertically to have the same
height as the constraining segment, but the result is a smoother looking band.

If `y` is a matrix, then banding is performed on (`t`,`y[:,i]`) for each column i.  The 
input `tmin_constraint` may be a vector with different entries for each column of `y`, or
can be a scalar value applied to all columns.  When `y` is a matrix, the outputs `band`,
`height`, and `tmin` are each vectors, one entry for each column of `y`.

# Examples #
```julia-repl
julia> t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
julia> (band,height,tmin) = linearband(t,y,0.4);
julia> height
0.3776184236760187
julia> p=plot(t,y,seriescolor=:blue,markersize=2,seriestype=:scatter)
julia> (LB,UB) = band
julia> plot!(p,LB[:,1],LB[:,2],seriescolor=:red)
julia> plot!(p,UB[:,1],UB[:,2],seriescolor=:red)
```

# Extended help

The algorithm is described in
Emily K. Szusz, Allan R. Willms, 2010 "A linear time algorithm for near minimax
continuous piecewise linear representations of discrete data", SIAM J. Sci.
Comput. 32 (5), pp. 2584-2602, doi=10.1137/090769077.

This implementation has made a few adjustments to the algorithm compared to the 
description in that paper.  In particular, segments whose slopes are between their
neighbours' slopes are allowed to be swallowed.  Also, if a violation of the 
constraints is encountered, the algorithm focuses on splitting nearby segments to
see if the violation can be reversed.  If not the minimun height is considered to
be the height of the band whose splitting caused the violation.  The algorithm may
then continue to split less tall bands, down to `minheightfrac` of the minimum height.

# Copyright 2022 (Julia version) Allan R. Willms #

LinBound is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

Dr. Allan R. Willms,
Dept. of Mathematics and Statistics,
University of Guelph,
Guelph, ON N1G 2W1,
Canada
"""
function linearband(t,y,tmin_constraint;minheightfrac=0.8,variableheight=false)

	# Number of data points.
	n = length(t)
	if length(y) != n
		error("length of t and rows of y must be the same")
	end
	if tmin_constraint < 0.0
		error("tmin_constraint must be nonnegative")
	end
	# scale the data, both t and y, to be in [0,1]
	(tshift,tscale) = scale!(t)
	(yshift,yscale) = scale!(y)
	tmin_constraint /= tscale

	# Determine the minimum distance between data points and compare with tmin_constraint.
	if n < 3 || tmin_constraint <= minimum(t[2:n] .- t[1:n-1])
		# optimal solution is simply to connect the data
		return (([t y],[t y]),0.0,tmincalc0(t,y))
	end

	# maxseg is the maximum number of segments allowed in a solution
	maxseg = n - 1
	# allseg keeps track of the total number of elements in the segment array.  This is
	# in general larger than the total number of segments since each time a segment is 
	# split, the old one is kept in the array and two new ones are added to the end.
	allseg = 1
	########### Initialization of the first segment and pivot points ##########
	#
	# The segments are structures containing the absolute indices of the first 
	# and last data points within the segment, the slope of the band around the 
	# segment, the height of the segment, an array of 3 pivot point locations 
	# (relative to the segment) and an array of 3 pivottypes.

	# The first, middle, and last points are the initial pivots.  Their types are 
	# unknown=0.
	segmentarraysize = 2^ceil(Int64,log(n)/log(2))
	segment = Vector{Segment}(undef,segmentarraysize)
	segment[1] = Segment(1,n,true,0.0,0.0,MVector{3,Int64}(1,div(n+1,2),n),
						 MVector{3,Int64}(zeros(Int64,3)))
	# Process the first segment to minimize its height and determine its slope
	segment[1] = sandwich!(segment[1],t,y)

	# Initialize the array of intermediate solution structures.
	solnarraysize = 20
	soln = Vector{Soln}(undef,solnarraysize)
	soln[1] = Soln([1],[t[1],t[n]],[0.0,0.0],segment[1].height,Inf)

	# There are three Modes:
	# reduction : the current solution has tmin >= tmin_constraint, split the tallest
	#	  segment that is not marked as unsplittable, and has positive height, until no
	#	  such segments remain.  If recovery mode has been entered once and failed, then
	#	  the min height of the solution has been found.  In this case continue splitting
	#	  the next tallest segments that are below minheightfrac*soln.height, 
	#	  or would have split lengths (in terms of number of data points) long enough.???
	# recovery : the current solution has tmin < tmin_constraint, split segments near
	#     the violating segment to try to alter this situation, if unsuccessful, mark the
	#     violating segment as unsplittable, return to reduction mode
	# smoothing : go through all segments splitting those that are around corners.
	mode = reduction
	count = 0
	current = 1  # index into soln array of the current solution being worked on
	seg2split = 1
	recoverlevel = 0
	# recoverlevelmax is the maximum number of times we try splitting neighbouring
	# segments that are not helpful when working in recovery mode
	recoverlevelmax = 4
	recover = Vector{Tuple{Int64,Int64,Int64}}(undef,0)
	violatingsegment = 0
	while seg2split != 0
		println("   going to split segment $(seg2split) of $(length(soln[current].segsort)) whose height is $(segment[soln[current].segsort[seg2split]].height)")
		# allocate more segment space if needed
		if allseg > segmentarraysize - 2
			segmentarraysize += n
			resize!(segment,segmentarraysize)
		end
		# split the selected segment of the current solution
		segment[allseg .+ [1,2]] = splitseg(segment[soln[current].segsort[seg2split]])
		# sandwich both new segments
		for i=1:2
			segment[allseg + i] = sandwich!(segment[allseg + i],t,y)
		end
		allseg += 2
		# construct the segsort and tjoin arrays for the new solution
		newsegsort = [soln[current].segsort[1:seg2split-1]; allseg-1; allseg; 
					  soln[current].segsort[seg2split+1:end]]
		newtjoin = [soln[current].tjoin[1:seg2split]; NaN; soln[current].tjoin[seg2split+1:end]]
		newyjoin = [soln[current].yjoin[1:seg2split]; NaN; soln[current].yjoin[seg2split+1:end]]
		# determine which tjoin points must be re-computed.  At least seg2split and the
		# two following, but maybe more to the left.
		index = segmentstojoin(newtjoin,seg2split)
		println(" Going to compute tjoin for indices $(index[1]) to $(index[2])")
		println(" currently they are:")
		display(newtjoin[index[1]:min(length(newtjoin),index[2])])
		# compute the tjoin points for the indices from index[1] to index[2]
		joinsegs!(newtjoin,newyjoin,index,newsegsort,segment,t,y)
		# compute tmin and the tjoinsw points
		(tmin, tminloc) = computetmin(newtjoin,newsegsort,segment)
		if tmin < tmin_constraint
			@show tminloc
		end
		# find the height of the new solution 
		height = computeheight(newtjoin,newsegsort,segment)

		if tmin >= tmin_constraint
			# This is the best solution so far, replace the solution array with just this
			# solution.
			soln[1] = Soln(newsegsort,newtjoin,newyjoin,height,tmin)
			count += 1
			println("solution $(count)")
			if mode == recovery
				println("recovery mode worked! at recoverlevel = $(recoverlevel) current = $(current)")
				#(LB,UB) = calculateband(soln[1],segment,t,y,variableheight,minheightfrac)
				#(err,val,ind) = testband(tmin_constraint,LB)
				#if err
				#	@show soln[1].tjoin[ind-2:ind+4]
				#	error("error val is $val but tmin is $tmin")
				#end
			end
			current = 1
		else
			# This solution is invalid but may yet be recovered, record to solution array.
			if current == solnarraysize  # make bigger if necessary
				solnarraysize += 20
				resize!(soln,solnarraysize)
			end
			current = current + 1
			soln[current] = Soln(newsegsort,newtjoin,newyjoin,height,tmin)
		end
		# Determine if we need to change modes and next segment to split.
		if mode != smoothing
			if tmin >= tmin_constraint
				if mode == recovery
					splice!(recover,1:length(recover))  # empty the recover array
					recoverlevel = 0
					mode = reduction
				end
				seg2split = findnextreductionsplit(soln[1],segment,minheightfrac)
				if seg2split == 0 || length(soln[current].segsort) == maxseg || height == 0.0
					println("Entering smoothing mode")
					mode = smoothing   # move to smoothing mode
					seg2split = findnextsmoothingsplit(soln[1],segment,2)
				end
			else
				if mode == reduction
					violatingsegment = seg2split  # record the violating segment
					mode = recovery
					println("      violator: segment = $(violatingsegment)")
				end
				recoverlevel += 1
				println("         recover level now $(recoverlevel)")
				if recoverlevel <= recoverlevelmax  # find next split
					segs = findnextrecoverysplit(soln[current],segment,tminloc)
					# If both entries of segs are not 0, then the first will be split
					# now and the second added to the recover array as a possible train.
					seg2split = segs[1]
					if segs[2] != 0  
						push!(recover,(recoverlevel,current,segs[2]))
					end
				end
				# If recoverlevel is above max or both entries of segs are 0, then 
				# nothing more can be done on this train.
				if recoverlevel > recoverlevelmax || segs[1] == 0 # give up on this train, pop another 
					if length(recover) > 0
						(recoverlevel,current,seg2split) = popfirst!(recover)
					else  # no recovery trains left, record violating segment as unsplittable
						println("      recovery failed, marking segment $(violatingsegment) as unsplittable")
						mode = reduction
						current = 1
						recoverlevel = 0
						# reduce the max recover level since we now have the constraining
						# segment and we don't wish to spend too much time on trying to
						# recover from further violating splits
						recoverlevelmax = 1
						segment[soln[1].segsort[violatingsegment]] = 
								makeunsplittable(segment[soln[1].segsort[violatingsegment]])
						seg2split = findnextreductionsplit(soln[1],segment,minheightfrac)
						if seg2split == 0 || length(soln[1].segsort) == maxseg || soln[1].height == 0.0
							println("Entering smoothing mode from recovery")
							println("There are $(length(soln[1].segsort)) segments in the solution")
							tjoincheck(soln[1].tjoin)
							mode = smoothing   # move to smoothing mode
							seg2split = findnextsmoothingsplit(soln[1],segment,2)
						end
					end
				end
			end
		else   # mode == smoothing 
			tjoincheck(soln[1].tjoin)
			# if the split was successful just find the next to split
			##@show [segment[soln[1].segsort[j]].slope for j = seg2split-1:seg2split+2]
			##@show soln[1].tjoin[seg2split-1:seg2split+2]
			if tmin >= tmin_constraint
				seg2split = findnextsmoothingsplit(soln[1],segment,seg2split)
			else  # unsuccessful split, so ignore new solution and find the next segment to split
				current = 1
				seg2split = findnextsmoothingsplit(soln[1],segment,seg2split+1)
			end
		end
	end
	# Calculate (t,y) values for the lower and upper boundaries of the band
	println("tjoin at end of main loop:")
	@show soln[1].tjoin[25:29]
	println("t starts of segments 24 to 29")
	@show [t[segment[soln[1].segsort[j]].first] for j in 24:29]
	tt = soln[1].tjoin[.!isnan.(soln[1].tjoin)]
	println("minimum delta tjoin is: $(minimum(tt[2:end] .- tt[1:end-1]))")
	println(" tdiff is $(t[56]-t[49])")
	(LB,UB) = calculateband(soln[1],segment,t,y,tmin_constraint,variableheight,minheightfrac)

	println("testing lower band")
	testband(tmin_constraint,LB,soln[1],segment)
	println("testing upper band")
	testband(tmin_constraint,UB,soln[1],segment)
	# rescale back to original units
	@views begin
		LB[:,1] .= tscale.*LB[:,1] .+ tshift
		LB[:,2] .= yscale.*LB[:,2] .+ yshift
		UB[:,1] .= tscale.*UB[:,1] .+ tshift
		UB[:,2] .= yscale.*UB[:,2] .+ yshift
	end
	return ((LB,UB),soln[1].height*yscale,soln[1].tmin*tscale)
end

function linearband(t, y::Matrix{T} where T<:Number, tmin_constraint; minheightfrac=0.8, variableheight=false)
	n = size(y,2)
	if length(tmin_constraint) == 1
		tmin_constraint = fill(tmin_constraint[1], n)
	end
	band = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, n)
	height = Vector{Float64}(undef,n)
	tmin = Vector{Float64}(undef,n)
	for i in 1:n
		(band[i],height[i],tmin[i]) = linearband(t, view(y,:,i), tmin_constraint[i]; 
											   minheightfrac=minheightfrac, variableheight=variableheight)
	end
	return (band, height, tmin)
end

function makeunsplittable(seg)
	# replace the segment seg with an equivlant one that is marked unsplittable
	return Segment(seg.first,seg.last,false,seg.slope,seg.height,seg.pivot,seg.pivottype)
end

function calculateband(soln,segment,t,y,tmin_constraint,variableheight,minheightfrac)
	# calculate the lower bound, LB, and upper bound, UB, of a band defined by this solution.
	# The first column of both LB and UB are the t values, and the second columns are the y values.
	# In the variable height case the t values will be different.
	if variableheight
		(LBtjoin,LByjoin,UBtjoin,UByjoin) = calculatevariabletjoin(soln,segment,t,y,tmin_constraint,minheightfrac)
		LB = calculatebound(-1,LBtjoin,LByjoin,0.0)
		UB = calculatebound(1,UBtjoin,UByjoin,0.0)
	else
		LB = calculatebound(-1,soln.tjoin,soln.yjoin,soln.height)
		UB = calculatebound(1,soln.tjoin,soln.yjoin,soln.height)
	end
	return (LB,UB)
end

function calculatevariabletjoin(soln,segment,t,y,tmin_constraint,minheightfrac)
	# calculate the tjoin and yjoin values for the variable height band
	println("Recomputing tjoin and yjoin LB and UB")
	n = length(soln.tjoin)
	println(" len tjoin = $(n),  len segsort = $(length(soln.segsort))")
	LBtjoin = Vector{Float64}(undef,n)
	LByjoin = Vector{Float64}(undef,n)
	UBtjoin = Vector{Float64}(undef,n)
	UByjoin = Vector{Float64}(undef,n)
	# Recalculate all tjoin points for the variable height case.  We do this one index at 
	# at time from left to right because it is possible that some segments heights are 
	# updated by one bound calculation and this will affect the calculation for the other bound.
	jlb = jub = i = 1
	while true
		jlb = joinseg_bound!(LBtjoin,LByjoin,i,-1,soln.segsort,segment,t,y,tmin_constraint)
		if (i = min(jlb,jub)) > n
			break
		end
		jub = joinseg_bound!(UBtjoin,UByjoin,i,1,soln.segsort,segment,t,y,tmin_constraint)
		if (i = min(jlb,jub)) > n
			break
		end
	end
	return (LBtjoin,LByjoin,UBtjoin,UByjoin)
end

function calculatebound(boundtype,tjoin,yjoin,height)
	# calculate the lower or upper bound from the tjoin and yjoin values
	n = sum(.!isnan.(tjoin))
	bnd = Matrix{Float64}(undef,n,2)
	j = 0
	for i = 1:n
		j = findnext(!isnan, tjoin, j+1)
		bnd[i,:] = [tjoin[j], yjoin[j] + boundtype*height*0.5]
	end
	return bnd
end

function tmincalc0(t,y)
	# compute tmin in the case that the optimal solution is simply to connect all data
	# points in a zero height band.
	slopecalc(i) = (y[i] - y[i-1])/(t[i] - t[i-1])
	slope = [0.0, slopecalc(2), slopecalc(3)]
	tmin = Inf
	for i=4:length(t)
		slope[1:3] = [slope[2:3];slopecalc(i)]
		if xor(slope[1] < slope[2], slope[2] < slope[3])
			len = t[i-1] - t[i-2]
			if len < tmin
				tmin = len
			end
		end
	end
	return tmin
end

function testband(tmin_constraint,B,soln,segment)
	err = false
	errval = 0.0
	errind = 0
	n = size(B,1)
	slope = Vector{Float64}(undef,n-1)
	for i=1:n-1
		slope[i] = (B[i+1,2] - B[i,2])/(B[i+1,1] - B[i,1])
	end
	#display([soln.tjoin[.!isnan.(soln.tjoin)] B[:,1]])
	#nseg = length(soln.segsort)
	#ind = findall(!isnan, soln.tjoin)
	#if ind[end] > nseg
	#	pop!(ind)
	#end
	#display([[segment[soln.segsort[j]].slope for j in ind] slope])
	for i=2:n-2
		val = B[i+1,1] - B[i,1]
		if val - tmin_constraint < -tmin_constraint*1.0e-14 && 
			!((slope[i-1] < slope[i] < slope[i+1]) ||
			 (slope[i-1] > slope[i] > slope[i+1]))
			println("bound is invalid at i = $(i), t = [$(B[i,1]),$(B[i+1,1])]")
			println("B[$(i-1):$(i+1),2] =")
			@show B[max(1,i-2):min(n,i+2),:]
			@show slope[max(1,i-2):min(n,i+1)]
			errval = val
			errind = i
			err = true
			break
		end
	end
	return (err,errval,errind)
end
		
function tjoincheck(tjoin)
	left = findfirst(!isnan,tjoin)
	right = findnext(!isnan,tjoin,left+1)
	while right != nothing
		val = tjoin[right] - tjoin[left]
		if val < 0.0
			println("  ************** tjoin not valid at $(left) $(right)")
			println("     tjoin[$(left)] = $(tjoin[left]), tjoin[$(right)] = $(tjoin[right])")
			error("tjoincheck failed")
		end
		left = right
		right = findnext(!isnan,tjoin,left+1)
	end
	return nothing
end

function scale!(x)
	# shift and scale the elements of x so that they are in [0,1].  
	#    shift = minimum(x)
	#    scale = maximum(x) - minimum(x)
	#    updated_x = (x - shift)/scale
	# Return the scale factor and the shift.
	shift = minimum(x)
	scale = maximum(x) - shift
	x .= (x .- shift)./scale
	return (shift,scale)
end

end
