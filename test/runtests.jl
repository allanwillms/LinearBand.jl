using LinearBand
using Test

@testset "LinearBand.jl" begin
    # Write your tests here.
	function tmintest(tmin_constraint,B)
		n = size(B,1)
		slope = Vector{Float64}(undef,n-1)
		for i=1:n-1
			slope[i] = (B[i+1,2] - B[i,2])/(B[i+1,1] - B[i,1])
		end
		for i=2:n-2
			val = B[i+1,1] - B[i,1]
			if val - tmin_constraint < -tmin_constraint*1.0e-14 && 
				!((slope[i-1] < slope[i] < slope[i+1]) ||
				 (slope[i-1] > slope[i] > slope[i+1]))
				return false
			end
		end
		return true
	end
	@testset "input consistency tests" begin
		t = range(0.0,10.0,10); y = ones(10);
		@test_throws ArgumentError (band,height,tmin) = LinearBand.linearband(t,y,-0.9);
		y = ones(5);
		@test_throws DimensionMismatch (band,height,tmin) = LinearBand.linearband(t,y,0.9);
	end
	@testset "tmin_constraint smaller than delta t" begin
		t = range(0.0,10.0,11); y = [1.0,2,4,3,6,5,3,4,1,8,2];
		(band,height,tmin) = LinearBand.linearband(t,y,0.9);
		@test height == 0.0
		@test tmin == 1.0
		@test all(band[1][:,1] .- t .== 0.0)
		@test all(band[2][:,1] .- t .== 0.0)
		@test all(band[1][:,2] .- y .== 0.0)
		@test all(band[2][:,2] .- y .== 0.0)
	end
	@testset "Constant height band" begin
		t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
		tmin_constraint = 0.4;
		(band,height,tmin) = LinearBand.linearband(t,y,tmin_constraint);
		@test height ≈ 0.377618 atol = 1.0e-6
		@test tmin ≈ 1.362962 atol = 1.0e-6
		@test band[1][1,1] == t[1]
		@test band[2][1,1] == t[1]
		@test band[1][end,1] == t[end]
		@test band[2][end,1] == t[end]
		@test all(band[1][:,1] .== band[2][:,1])
		@test all(band[1][:,2] .<= band[2][:,2])
		@test all(band[1][2:end,1] .- band[1][1:end-1,1] .> 0.0)
		@test band[1][5,2] ≈ -1.011158 atol = 1.0e-6
		@test band[2][7,2] ≈ -0.381188 atol = 1.0e-6
		@test tmintest(tmin_constraint,band[1])
		@test tmintest(tmin_constraint,band[2])
	end
	@testset "Variable height band" begin
		t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
		tmin_constraint = 0.3;
		(band,height,tmin) = LinearBand.linearband(t,y,tmin_constraint,variableheight=true);
		@test height ≈ 0.178072 atol = 1.0e-6
		@test tmin ≈ 0.341476 atol = 1.0e-6
		@test band[1][1,1] == t[1]
		@test band[2][1,1] == t[1]
		@test band[1][end,1] == t[end]
		@test band[2][end,1] == t[end]
		@test band[1][15,1] - band[2][15,1] ≈ 0.026049 atol = 1.0e-6
		@test band[1][3,2] ≈ 0.890467 atol = 1.0e-6
		@test band[2][21,2] ≈ 0.694347 atol = 1.0e-6
		@test tmintest(tmin_constraint,band[1])
		@test tmintest(tmin_constraint,band[2])
	end
	@testset "Matrix input" begin
		t = range(0.0,10.0,50); y = zeros(50,2); y[:,1] = @. sin(t) + 0.2*cos(9*t);
		y[:,2] = @. exp(-0.05*t)*cos(t) + 0.15*sin(13*t);
		(band,height,tmin) = LinearBand.linearband(t,y,0.35);
		@test length(band) == 2
		@test length(height) == 2
		@test length(tmin) == 2
		@test height[1] ≈ 0.347954 atol = 1.0e-6
		@test height[2] ≈ 0.271780 atol = 1.0e-6
		@test tmin[1] ≈ 0.988274 atol = 1.0e-6
		@test tmin[2] ≈ 0.381371 atol = 1.0e-6
		@test band[1][1][11,2] ≈ 0.830181 atol = 1.0e-6
		@test band[2][2][8,1] ≈ 4.388650 atol = 1.0e-6
	end
end
