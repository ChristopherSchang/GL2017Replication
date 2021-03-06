using Test, GL2017Replication, MAT

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
    initilize!(gl)
    @test sum(gl.pr) ≈ 1.0 
end

@testset "EGM" begin
    gl = ModelGL()
    initilize!(gl)
    EGM!(gl)
    @test minimum(gl.b_pol)>=gl.b_grid[1]
    @test maximum(gl.b_pol)<=gl.b_grid[end]
end

@testset "compute distribution" begin
    gl = ModelGL()
    initilize!(gl)
    EGM!(gl)
    compute_distribution!(gl)
    @test sum(gl.JD) ≈ 1
end

@testset "aggregate" begin
    gl = ModelGL()
    initilize!(gl)
    EGM!(gl)
    compute_distribution!(gl)
    aggregate!(gl)
    @test gl.D_4Y_actual < gl.B_4Y_actual
end

cd(@__DIR__)
steady_mat = matopen("steady.mat")
@testset "calibrate" begin
    gl = ModelGL()
    load_parameters_iss!(gl)
    calibrate!(gl)
    @testset "initial" begin
        Y1_mat =  read(steady_mat, "Y1")
        cpol1_mat = read(steady_mat, "c_pol1")
        @test abs(gl.Y_actual - Y1_mat) < 0.001
        @test abs(maximum(gl.c_pol) - maximum(cpol1_mat)) < 0.001
    end 
 
    gl_tss = calibrate_terminal(gl,load_params = true)

    @testset "calibrate_terminal" begin  
        Y2_mat =  read(steady_mat, "Y2")
        cpol2_mat = read(steady_mat, "c_pol2")
        @test abs(gl_tss.Y_actual - Y2_mat) < 0.001
        @test abs(maximum(gl_tss.c_pol) - maximum(cpol2_mat)) < 0.001
    end

    #Takes very long to run 
    #Tgl = TransGL()
    #gl_trans = transition!(gl,gl_tss,Tgl)
    #
    #@testset "transition" begin   
    #    @test abs(Tgl.D_4Y_t[1] - gl.D_4Y_actual) < 0.001
    #    @test abs(Tgl.D_4Y_t[end] - gl_tss.D_4Y_actual) < 0.001
    #end
     
end
 