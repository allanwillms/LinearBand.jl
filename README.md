# LinearBand

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

---

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

[![Build Status](https://github.com/allanwillms/LinearBand.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/allanwillms/LinearBand.jl/actions/workflows/CI.yml?query=branch%3Amain)
