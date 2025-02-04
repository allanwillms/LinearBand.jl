using LinearBand
using Test

@testset "LinearBand.jl" begin
    # Write your tests here.
	@testset "Constant height band" begin
		t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
		(band,height,tmin) = LinearBand.linearband(t,y,0.4)
		@test height ≈ 0.377618 atol = 1.0e-6
		@test tmin ≈ 1.362962 atol = 1.0e-6
		@test band[1][1,1] == t[1]
		@test band[2][1,1] == t[1]
		@test band[1][end,1] == t[end]
		@test band[2][end,1] == t[end]
		@test all(band[1][:,1] .== band[2][:,1])
		@test all(band[1][:,2] .<= band[2][:,2])
		@test band[1][5,2] ≈ -1.011158 atol = 1.0e-6
		@test band[2][7,2] ≈ -0.381188 atol = 1.0e-6
	end
	@testset "Variable height band" begin
		t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
		(band,height,tmin) = LinearBand.linearband(t,y,0.3,variableheight=true)
		@test height ≈ 0.106988 atol = 1.0e-6
		@test tmin ≈ 0.387886 atol = 1.0e-6
		@test band[1][1,1] == t[1]
		@test band[2][1,1] == t[1]
		@test band[1][end,1] == t[end]
		@test band[2][end,1] == t[end]
		@test band[1][15,1] - band[2][15,1] ≈ -0.005886 atol = 1.0e-6
		@test band[1][3,2] ≈ 0.864069 atol = 1.0e-6
		@test band[2][25,2] ≈ 0.693300 atol = 1.0e-6
	end
	@testset "Matrix input" begin
		t = range(0.0,10.0,50); y = zeros(50,2); y[:,1] = @. sin(t) + 0.2*cos(9*t);
		y[:,2] = @. exp(-0.05*t)*cos(t) + 0.15*sin(13*t)
		(band,height,tmin) = LinearBand.linearband(t,y,0.35)
		@test length(band) == 2
		@test length(height) == 2
		@test length(tmin) == 2
		@test height[1] ≈ 0.347954 atol = 1.0e-6
		@test height[2] ≈ 0.272982 atol = 1.0e-6
		@test tmin[1] ≈ 0.988273 atol = 1.0e-6
		@test tmin[2] ≈ 0.666448 atol = 1.0e-6
		@test band[1][1][11,2] ≈ 0.830181 atol = 1.0e-6
		@test band[2][2][8,1] ≈ 5.926001 atol = 1.0e-6
	end
end
