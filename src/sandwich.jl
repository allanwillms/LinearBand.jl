function sandwich!(segment, t, y)
	# Apply the sandwich (Remez) algorithm to the data segment.
	# It alters pivot field of segment structure, but the expectation
	# is that you are going to overwrite the segment anyway.
	segind = segment.first:segment.last
	# if there only two points in the segment, then height is zero and res is zero
	if length(segind) < 3
		slope = (y[segment.last] - y[segment.first])/(t[segment.last] - t[segment.first])
		return Segment(segment.first,segment.last,segment.splittable,slope,0.0,[1,2,2],[0,0,0])
	end
	tseg = t[segind]
	yseg = y[segind]
	(exchangemade,res,E,slope) = remez1step!(segment.pivot, tseg, yseg)
	n = 1
	while exchangemade
		n = n+1
		(exchangemade,res,E,slope) = remez1step!(segment.pivot, tseg, yseg)
		if n == 20
			error("took 20 iterations of remez")
		end
	end
	return Segment(segment.first,segment.last,segment.splittable,slope,2.0*E,segment.pivot,Int64.(sign.(res)))
end

function remez1step!(pivot, t, y)
	# Apply one iteration of the Remez algorithm to the data segment.
	# Pivot locations are altered.  Return residues at the new pivots and the value of E.
	#
	# (t,y) are the points in one segment, pivot is an array of increasing indices
	# into (t,y) relative to this segment (so 1 is first index of t and y)
	# determine the best fit line through the three pivot points
	coeff = linefit(t[pivot],y[pivot])
	E = abs(coeff[3])
	# compute the residuals at all points of the segment
	res = y .- (coeff[1] .+ coeff[2].*t)
	# exchange pivots
	exchangemade = exchangepivots!(pivot,res,E,y)
	return (exchangemade,res[pivot],E,coeff[2])
end

function linefit(t,y)
	# solve the three equations  y_i = c0 + c1 * t_i + (-1)^i * E,  i=1,2,3
	# for c0, c1, and E.
	c1 = (y[3] - y[1])/(t[3] - t[1])
	temp1 = y[2] - t[2]*c1
	temp2 = y[1] - t[1]*c1
	return [(temp1 + temp2)*0.5, c1, (temp1 - temp2)*0.5]
end

function exchangepivots!(pivot,res,E,y)
	# The points y are just the points within one segment.  They have corresponding residues,
	# res.  pivot is an array of pivot locations relative to this segment, in increasing
	# order.  On input, the residues at all 3 pivot points will be the same (up to round
	# off error) absolute value |E| but will be alternating in sign.  Find the three
	# points with the largest magnitude residues.  For each, if its residue exceeds |E| in
	# magnitude (that is, exceeds by more than round off error) and if the magnitude of
	# its residue exceeds that of the residue of the nearest pivot point that has the same
	# sign residue, then delete that pivot point and make the new point a pivot point.  In
	# the event that the new point is to the right of all pivots and the nearest pivot has
	# opposite sign residue, then if the furthest left pivot point has smaller magnitude
	# residue then delete that pivot and make the new point a pivot.  (similarly if the
	# new point is to the left of all pivot points).  Ensure that pivots at end are sorted
	# in increasing order.  Returns true if a pivot exchange is made.
	#

	roundoff = 2.0*maximum(abs.(E .- abs.(res[pivot])))  # each of |res[pivot[i]]| should be E
	I = partialsortperm(res, 1:3; by=abs, rev=true)  # finds the locations of three largest |res|
	if sort(I) == sort(pivot)
		return false
	end
	exchangemade = false
	for ind in I
		absresind = abs(res[ind])
		if absresind > E + roundoff
			signresind = sign(res[ind])
			# determine which current pivot to compare with ind
			if ind < pivot[1]   # left of leftmost pivot
				if signresind == sign(res[pivot[1]])     # same sign res as pivot[1]
					compare = 1
				else    # opposite sign res as pivot[1]
					compare = 3
				end
			elseif ind > pivot[3]   # right of rightmost pivot
				if signresind == sign(res[pivot[3]])     # same sign res as pivot[3]
					compare = 3
				else    # opposite sign res as pivot[3]
					compare = 1
				end
			else # between left and right pivots
				if signresind == sign(res[pivot[2]])    # same sign as pivot[2]
					compare = 2
				else      # opposite sign as middle pivot, so same sign as left and right
					if ind < pivot[2]
						compare = 1
					else 
						compare = 3
					end
				end
			end
			# Since on input all pivots have same residue, the first exchange will always
			# satisfy the below test, but the test is necessary since the second and third
			# largest points may be bigger than |E| but not bigger than largest and so we do
			# want to exchange them with the largest.
			if absresind > abs(res[pivot[compare]]) + roundoff 
				exchangemade = true
				pivot[compare] = ind
				sort!(pivot)  # need to sort to maintain increasing order
			end
		end
	end
	return exchangemade
end
