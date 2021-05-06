using Test, GL2017Replication

@testset "Testing the testing" begin
    @test hello("Julia") == "Hello, Julia"
    @test domath(2.0) ≈ 7.0 
end

@testset "Structures" begin
    gl = ModelGL()
    @test typeof(gl.β) == Float64 
end

@testset "initilize" begin
    gl = ModelGL() 
    prior_pr = copy(gl.pr)
    initilize!(gl)
    @test sum(gl.pr) ≈ 1.0 
end