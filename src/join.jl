function joinsegs!(tjoin,yjoin,indices,segsort,segment,t,y)
	# call joinseg! (boundtype=0 (midline)) for each i in indices[1]:indices[2]
	# After each call, test to see if the segment to the left must be swallowed.
	# If so, we will need to recompute tjoin.
	for i in indices[1]:indices[2]
		joinseg!(tjoin,yjoin,i,0,segsort,segment,t,y)
	end
end

function joinseg!(tjoin,yjoin,i,boundtype,segsort,segment,t,y)
	# Update tjoin and yjoin values for the location i.  tjoin[i+1] is the
	# intersection of the lower bound, midline, or upper bound (for boundtype = 
	# -1, 0, or 1, resp.) of the neighbouring segment to the left and segment i+1 (thus
	# the left edge of segment i+1).  Usually the left neighbour is segment i, but segment
	# i may have been swallowed legally (its slope is between the slopes of its neighbours
	# and tjoin[i+1] was originally computed as being smaller than tjoin[i]), in which
	# case the left neighbour is one further to the left (i-1), unless it was swallowed,
	# etc.  If segment i is swallowed then we set tjoin[i] to NaN.  The first and last
	# non-NaN values of tjoin are t[1] and t[end], respectively.
	nseg = length(segsort)
	pt = Vector{Tuple{Float64,Float64}}(undef,2)
	left = findprev(!isnan, tjoin, i-1)
	finished = false
	while !finished
		# find the intersection of segment left and segment i
		if isnothing(left)  # i is the leftmost non-swallowed segment
			tjoin[i] = t[1]
			yjoin[i] = calculate_yval(t[1],getpoint(segment[segsort[i]],boundtype,1,t,y),
									  segment[segsort[i]].slope)
			finished = true
		else
			if i == nseg + 1
				tjoin[i] = t[end]
				yjoin[i] = calculate_yval(t[end],getpoint(segment[segsort[left]],3,boundtype,t,y),
										  segment[segsort[left]].slope)
			else
				(tjoin[i],yjoin[i]) = calc_seg_intersection!(pt,segment[segsort[[left,i]]],boundtype,t,y)
				#println("      yjoin[$(i)] is $(yjoin[i]),  left = $(left), i = $(i)")
				#display(segment[segsort[left]])
				#display(segment[segsort[i]])
				#println("")
			end
			if tjoin[i] - tjoin[left] <= 0.0  && (left == 1 || left == nseg ||
												  middleslope(segment,segsort,left-1:left+1))
				# legal swallow of segment left
				println("        legally SWALLOWING segment $(left) of $(nseg)")
				println("           slopes are: ")
				@show [segment[segsort[j]].slope for j in max(1,left-1):min(nseg,left+1)]
				println("           tjoins are: ")
				@show tjoin[[left,i]]
				tjoin[left] = NaN
				left = findprev(!isnan, tjoin, left-1)
			else
				println("Computed boundtype $(boundtype) tjoin[$(i)] as $(tjoin[i]), yjoin as $(yjoin[i])")
				println("   heights of segments were $(segment[segsort[left]].height), $(segment[segsort[min(nseg,i)]].height)")
				println("   delta is $(tjoin[i] - tjoin[left]), slopes are $([segment[segsort[j]].slope for j in max(1,left-1):min(length(segsort),left+1)])")
				if i != left+1 && i<=length(segsort)
					println("     slope of segment i=$(i) is $(segment[segsort[i]].slope)")
				end
				finished = true
			end
		end
	end
	return nothing
end

function joinseg_bound!(tjoin,yjoin,i,boundtype,segsort,segment,t,y,tmin_constraint)
	# Join the lower or upper bound for segment i and its left unswallowed neighbour.
	# This function is called after a valid set of tjoin points for the midline has been
	# found defining a valid constant height band.  This function is then called 
	# if the user wanted a variable height band.  In this case, it is possible that 
	# other segments could be legally swallowed (slope between neighbouring slopes) even
	# though the midline was not swallowed.  Such is taken care of by calling tjoin!.  It
	# is also possible that tmin_constraint violoations or non-legal swallowings may occur.  In this case we expand the
	# heights of the segment and its two neighbours (starting with the shortest and when
	# it reaches the next height, both are increased, etc.) until the tmin_constraint
	# condition is satisfied.
	# Expanding the heights of the segments will work to remove the constraint violation because once all three heights are the same the
	# tjoin points will be the same as the midline, and it was legal.
	# If the height of the band to the left or the middle is increased, then we will need
	# to recompute that tjoin (tjoins are computed left to right).
	# Return value is the next tjoin value to be computed.  This will either be i+1, or,
	# if any segment from i or lower has had its height increased, then it will be the
	# minimum index for such segments.
	nseg = length(segsort)
	joinseg!(tjoin,yjoin,i,boundtype,segsort,segment,t,y)
	middle = findprev(!isnan, tjoin, i-1)
	delta_t = 0.0
	if !isnothing(middle)
		left = findprev(!isnan, tjoin, middle-1)
		if !isnothing(left) && i < nseg && !middleslope(segment,segsort,[middle-1,middle,middle+1])
			delta_t = tjoin[i] - tjoin[middle] - tmin_constraint
		end
	else 
		left = nothing
	end
	println(" i=$(i), middle=$(middle), left=$(left), delta_t = $(delta_t)")
	if !isnothing(middle)
		@show segment[segsort[middle]]
		@show t[segment[segsort[middle]].first]
	end
	if i <= nseg
		@show segment[segsort[i]]
		@show t[segment[segsort[i]].first]
	end
	println(" computed tjoin, yjoin as: $(tjoin[i]), $(yjoin[i])")
	if !isnothing(middle) && !isnothing(left) && i < nseg && delta_t < 0.0 #-tmin_constraint*1.0e-14
		# record heights
		ind = [left,middle,i]
		heights_orig = [segment[segsort[j]].height for j in ind]
		@show heights_orig
		p = sortperm(heights_orig)
		heights_sorted = heights_orig[p]
		@show heights_sorted
		slopes = [segment[segsort[j]].slope for j in ind]
		@show slopes
		# hfun is the change in height of each segment, j = 1,2, or 3.
		hfun(j,h) = max(heights_orig[j],h) - heights_orig[j]
		# joinfun will be the possible new tjoin point, j should be 2 or 3
		joinfun(j,h) = tjoin[ind[j]] + boundtype*0.5*(hfun(j-1,h) - hfun(j,h))/(slopes[j] - slopes[j-1])
		# delta is the time difference in join points for various heights
		# By construction, delta(heights_sorted[1]) = delta_t.
		# delta is a piecewise linear function of h
		delta(h) = joinfun(3,h) - joinfun(2,h) - tmin_constraint
		val2 = delta(heights_sorted[2])
		@show delta_t
		@show val2
		if val2 >= 0.0
			# we only need to raise the height of the lowest segment part way to
			# heights_sorted[2]
			# @show heights_sorted
			height = (heights_sorted[1]*val2 - heights_sorted[2]*delta_t)/(val2 - delta_t)
			println("updating height for segment $(ind[p[1]]) from $(segment[segsort[ind[p[1]]]].height) to $(height)")
			newheight = (heights_sorted[1]*val2 - heights_sorted[2]*delta_t)/(val2 - delta_t)
			# increase height by 1.0e-6 the distance between two lowest heights to avoid
			# round off errors causing the height to be bit too small
			newheight = min(heights_sorted[2],newheight*(1.0 + 1.0e-6*(heights_sorted[2] - heights_sorted[1])))
			updateheight!(segment,segsort[ind[p[1]]],newheight)
			nextcompute = ind[p[1]]
		else
			val3 = delta(heights_sorted[3])
			println("     val3 should be >0, is $(val3)")
			if val3 <= 0.0
				read(stdin,1)
			end
			# it must be that val3 > 0.  So height must be raised part way between heights_sorted[2]
			# and heights_sorted[3], so shortest two segments need updating.
			newheight = (heights_sorted[2]*val3 - heights_sorted[3]*val2)/(val3 - val2)
			# increase height by 1.0e-6 the distance between two largest heights to avoid
			# round off errors causing the height to be bit too small
			newheight = min(heights_sorted[3],newheight*(1.0 + 1.0e-6*(heights_sorted[3] - heights_sorted[2])))
			println("updating heights for segments $(ind[p[1]]) and $(ind[p[2]]) from $(segment[segsort[ind[p[1]]]].height) and $(segment[segsort[ind[p[2]]]].height) to $(newheight)")
			updateheight!(segment,segsort[ind[p[1]]],newheight)
			updateheight!(segment,segsort[ind[p[2]]],newheight)
			nextcompute = min(ind[p[1]],ind[p[2]])
		end
	else
		nextcompute = i+1
	end
	println("               ****************************")
	println("              nextcompute = $(nextcompute)")
	return nextcompute
end

function segmentstojoin(tjoin,loc)
	# tjoin values for loc to loc+2 must be recomputed, since loc and loc+1 are the new
	# segments resulting from the last split.  If tjoin[x] for x=loc-1,loc-2,... are
	# currently NaN (these segments were swallowed), then also recompute tjoin[x].
	# start is 1 right of the last nonNaN left of loc.
	# If tjoin[x] is NaN for x=loc+2,loc+3,... then also compute tjoin[x]
	# We also add one more tjoin to the right, even though that value will not be altered,
	# to ensure that that segment is tested for swallowing (which would then cause tjoin
	# to be altered).
	start = findlast(!isnan,@view(tjoin[1:loc-1]))
	if isnothing(start)
		start = 1
	else
		start += 1
	end
	n = length(tjoin)
	fin = findfirst(!isnan,@view(tjoin[loc+2:n]))
	if isnothing(fin)
		fin = n
	else
		loc += fin + 2
		fin = findfirst(!isnan,@view(tjoin[loc:n]))
		if isnothing(fin)
			fin = n
		else
			fin += loc-1
		end
	end
	return (start,fin)
end

function calc_seg_intersection!(pt,segments,boundtype,t,y)
	# calculate the intersection of the lower bounds, midlines, or upper bounds (for
	# boundtype = -1, 0, or 1, resp.) of two given segments (which should be ordered left
	# to right),
	# pt is an array of two tuples of 2 floats, used for storage
	for j=1:2
		pt[j] = getpoint(segments[j],boundtype,j,t,y)  # j=2 will be converted to 3 by getpoint
	end
	println("  slopes are:")
	display((segments[1].slope,segments[2].slope))
	println("   points are:")
	@show pt
	(tint, yint) = lineintersect(pt[1],segments[1].slope,pt[2],segments[2].slope)
	return (tint,yint)
end

function getpoint(seg,boundtype,k,t,y)
	# get a point on the lower bound, midline, or upper bound (for boundtype = -1, 0,
	# or 1, resp.) of the segment.  The t-value is at pivot[k] (k should be 1 or 3).  
	# The segment height value may be larger than the minimum band going through the
	# pivots.  This happens when joinset_bound! is called in the variable height case.  So
	# we compute the minimum band height from the first two pivot points.
	if k != 1
		k = 3
	end
	loc1 = seg.first - 1 + seg.pivot[k]
	loc2 = seg.first - 1 + seg.pivot[2]
	minbandheight = abs(y[loc2] - y[loc1] - seg.slope*(t[loc2] - t[loc1]))
	return (t[loc1], y[loc1] + (boundtype*seg.height - seg.pivottype[k]*minbandheight)*0.5)
end

function calculate_yval(t,pt,slope)
	# Given a value of t and a line defined by a point and a slope, calculate the value of
	# the line at t.
	return pt[2] + slope*(t - pt[1])
end

function lineintersect(pt1,slope1,pt2,slope2)
	# find the intersection of two lines, each specified by a point (t,y) and a slope y'.
	# returns (t_int, y_int)  If the slopes are the same but the lines are not coincident,
	# then tint with be +-Inf; if the lines are coincident then tint is set to the average
	# of the two input t values and yint computed accordingly.
	# We do this in a numerically stable manner.

	# First we compute the difference between the segment 2 line at t1 compared to y1
	differ1 = calculate_yval(pt1[1],pt2,slope2) - pt1[2]   
	# Now compute the difference between y2 and the segment 1 line at t2
	differ2 = pt2[2] - calculate_yval(pt2[1],pt1,slope1)
	if abs(differ1) < 1.0e-10 && abs(differ2) < 1.0e-10   # essentially coincident lines: average t values
		tint = 0.5*(pt1[1] + pt2[1])
	else
		rat = differ1/differ2
		if  abs(1.0 - rat) < 1.0e-10  # nearly parallel lines that do not cross between t1 and t2
			if abs(differ1) >= abs(differ2)
				tint = Inf*sign(pt2[1] - pt1[1])
			else
				tint = Inf*sign(pt1[1] - pt2[1])
			end
		else
			tint = pt1[1]/(1.0 - rat) + pt2[1]/(1.0 - 1.0/rat)
		end
	end
	yint = pt1[2] + slope1*(tint - pt1[1])
	return (tint, yint)
end

function computetmin(tjoin,segsort,segment)
	# Returns (tmin, tminloc), where tminloc is the location of the constraining segment.
	tmin = Inf 
	tminloc = 0
	nseg = length(segsort)
	left = findfirst(!isnan,tjoin)
	middle = findnext(!isnan,tjoin,left+1)
	right = findnext(!isnan,tjoin,middle+1)
	while right != nothing && right <= nseg
		seglen = tjoin[right] - tjoin[middle]
		if seglen < tmin && !middleslope(segment,segsort,[left,middle,right])
			tmin = seglen
			tminloc = middle
		end
		left = middle
		middle = right
		right = findnext(!isnan,tjoin,right+1)
	end
	return (tmin, tminloc)
end

function computeheight(tjoin,segsort,segment)
	# Returns height of the solution
	height = -1.0
	nseg = length(segsort)
	for i in 1:nseg
		if !isnan(tjoin[i]) && segment[segsort[i]].height > height
			height = segment[segsort[i]].height
		end
	end
	return height
end

function middleslope(segment,segsort,index)
	# return true if the slope of segement segsort[index[2]] lies strictly
	# between the slopes of segsort[index[1]] and segsort[index[3]]
	return (segment[segsort[index[1]]].slope < segment[segsort[index[2]]].slope <
			 segment[segsort[index[3]]].slope) ||
			(segment[segsort[index[1]]].slope > segment[segsort[index[2]]].slope >
			 segment[segsort[index[3]]].slope)
end

function updateheight!(segment,index,height)
	# update the height for segment index 
	segment[index] = Segment(segment[index].first,segment[index].last,segment[index].splittable,
							 segment[index].slope,height,segment[index].pivot,segment[index].pivottype)
	return nothing
end
