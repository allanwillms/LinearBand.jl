function splitseg(segment)
	# split the given segment at the middle pivot point
	# return an array of two segments
	# a new segment is unsplittable if it only has two points
	nleft = segment.pivot[2]
	breakpt = segment.first + nleft - 1
	leftpivots = [1,div(segment.pivot[2]+1,2),segment.pivot[2]]
	nright = segment.last - breakpt + 1
	rightpivots = [1,div(nright+1,2),nright]
	return [Segment(segment.first,breakpt,nleft>2,0.0,0.0,leftpivots,zeros(Int64,3)),
			Segment(breakpt,segment.last,nright>2,0.0,0.0,rightpivots,zeros(Int64,3))]
end

function findnextreductionsplit(soln,segment,minheightfrac)
	# find the tallest segment that is not marked as unsplittable
	# If their are no such segments or if the tallest such segment has height <= minheightfrac *
	# soln.height, then return 0
	f(i) = segment[soln.segsort[i]].splittable ? segment[soln.segsort[i]].height : -1.0
	(val,seg2split)= findmax(f, 1:length(soln.segsort))
	if val < minheightfrac*soln.height 
		seg2split = 0
	end
	return seg2split
end

function findnextrecoverysplit(soln,segment,viol)
	# determine the next segment to split when in recovery mode.  Segment viol violates
	# the constraint.  If a neighbouring segment is such that the slope discontinuity between 
	# viol and this neighbour is positive (negative) and the pivottype[2] for the neighbour is 
	# 1 (-1), then such a split will be helpful in that it will increase the length of segment 
	# viol.  In order to be splittable, the neighbour must not be marked as unsplittable.
	# Return a vector of the two neighbour segments ordered so that helpful neighbours are before 
	# unhelpful ones, and otherwise taller neigbours before shorter ones.  A value of 0 indicates
	# the neighbour is not splittable.
	increasing = [segment[soln.segsort[i+1]].slope - segment[soln.segsort[i]].slope > 0.0 
				  for i=viol-1:viol]
	ptype = [segment[soln.segsort[i]].pivottype[2] == 1 for i=[viol-1,viol+1]]
	helpful = @. !xor(increasing,ptype) # neighours whose split would help (1=left, 2=right)
	splittable = [segment[soln.segsort[i]].splittable for i=[viol-1,viol+1]]   
	if all(.!splittable)
		seg2split = [0,0]
	elseif all(splittable) 
		if !xor(helpful...)  # both helpful or both not helpful, choose tallest
			if segment[soln.segsort[viol-1]].height > segment[soln.segsort[viol+1]].height
				seg2split = [viol - 1,viol + 1]
			else
				seg2split = [viol + 1,viol - 1]
			end
		elseif helpful[1]
			seg2split = [viol - 1,viol + 1]
		else  # helpful[2]
			seg2split = [viol + 1,viol - 1]
		end
	elseif splittable[1]
		seg2split = [viol - 1,0]
	else # splittable[2]
		seg2split = [viol + 1,0]
	end
	return seg2split
end

function findnextsmoothingsplit(soln,segment,currentseg)
	# determine the segment starting at currentseg and moving right that is the next
	# segment for splitting in smoothing mode.  This means if must have a slope between that of its
	# neighbours and have the middle pivot be a lower pivot (if slopes increasing) or an upper
	# pivot (if slopes decreasing).  Also, the segment must have at least 3 points. 
	# Returns 0 if none are valid.
	while currentseg < length(soln.segsort)
		# determine if this segment is splittable
		increasing = [segment[soln.segsort[i+1]].slope - segment[soln.segsort[i]].slope > 0.0 
					  for i=currentseg-1:currentseg]
		if segment[soln.segsort[currentseg]].last - segment[soln.segsort[currentseg]].first > 1 &&
			((all(increasing) && segment[soln.segsort[currentseg]].pivottype[2] == -1) ||
			 (all(.!increasing) && segment[soln.segsort[currentseg]].pivottype[2] == 1))
			return currentseg
		else
			currentseg = currentseg + 1
		end
	end
	return 0
end
