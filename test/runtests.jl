using Test, JuliaDICE
using JuMP, Ipopt
using Pkg

res_cbopt = run_dice()

@testset "Test CB Optimal" begin
    @test res_cbopt.solved == true
    # ECO2
    @test isapprox(res_cbopt.ECO2[1], 42.9, rtol=0.1) 
    @test isapprox(res_cbopt.ECO2[2], 42.9, rtol=0.1) 
    @test isapprox(res_cbopt.ECO2[7], 37.1, rtol=0.1) 
    @test isapprox(res_cbopt.ECO2[17], 15.9, rtol=0.1)
    # MIU
    @test isapprox(res_cbopt.MIU[1], 0.05, rtol=0.01)
    @test isapprox(res_cbopt.MIU[3], 0.24, rtol=0.01)
    @test isapprox(res_cbopt.MIU[5], 0.31, rtol=0.01)
    @test isapprox(res_cbopt.MIU[7], 0.39, rtol=0.01)
    @test isapprox(res_cbopt.MIU[9], 0.46, rtol=0.01)
    @test isapprox(res_cbopt.MIU[17], 0.84, rtol=0.01)
end

res_crazy = run_dice(a2base = 0.01, bounds = Dict("MIU"=>("==",1.0), "TATM"=>("<=",15), "Y" =>(">=",[fill(floatmin(Float64),10);fill(0.1,71)]), "ECO2" =>("<=",10000)),optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0))

@testset "Test Crazy call" begin
    @test res_crazy.solved == false
end