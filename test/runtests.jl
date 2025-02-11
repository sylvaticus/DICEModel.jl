using Test, DICEModel
using JuMP, Ipopt


println("Testing DICEModel...")

res_cbopt    = run_dice_scenario("cbopt")
res_t2c      = run_dice_scenario("t2c")
res_t15c     = run_dice_scenario("t15c")
res_altdam   = run_dice_scenario("altdam")
res_parisext = run_dice_scenario("parisext")
res_base     = run_dice_scenario("base")
res_r5       = run_dice_scenario("r5")
res_r4       = run_dice_scenario("r4")
res_r3       = run_dice_scenario("r3")
res_r2       = run_dice_scenario("r2")
res_r1       = run_dice_scenario("r1")

@testset "Test CB Optimal" begin
    @test res_cbopt.solved == true
    # ECO2
    @test isapprox(res_cbopt.ECO2[1], 42.9, atol=0.051) 
    @test isapprox(res_cbopt.ECO2[2], 42.9, atol=0.051) 
    @test isapprox(res_cbopt.ECO2[7], 37.1, atol=0.051) 
    @test isapprox(res_cbopt.ECO2[17], 15.9, atol=0.051)
    # MIU
    @test isapprox(res_cbopt.MIU[1], 0.05, atol=0.0051)
    @test isapprox(res_cbopt.MIU[3], 0.24, atol=0.0051)
    @test isapprox(res_cbopt.MIU[5], 0.31, atol=0.0051)
    @test isapprox(res_cbopt.MIU[7], 0.39, atol=0.0051)
    @test isapprox(res_cbopt.MIU[9], 0.46, atol=0.0051)
    @test isapprox(res_cbopt.MIU[17], 0.84, atol=0.0051)
end

res_crazy = run_dice(a2base = 0.01, bounds = Dict("MIU"=>("==",1.0), "TATM"=>("<=",15), "Y" =>(">=",[fill(floatmin(Float64),10);fill(0.1,71)]), "ECO2" =>("<=",10000)),optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0))

@testset "Test Crazy call" begin
    @test res_crazy.solved == false
end

# Testing first and last values
# 5 vars, 10 scen --> 2 scen for each var

@testset "CO2 emissions with T<2°C and T<1.5°C scenarios" begin
    @test isapprox(res_t2c.ECO2[1], 42.9, atol=0.051)
    @test isapprox(res_t2c.ECO2[17], 1.2, atol=0.051)
    @test isapprox(res_t15c.ECO2[1], 42.9, atol=0.051)
    @test isapprox(res_t15c.ECO2[7], 5.7, atol=0.051)

end

@testset "CO2 concentration with alt damage and Paris extended scenarios" begin
    @test isapprox(res_altdam.ppm[1], 416.2, atol=0.1)
    @test isapprox(res_altdam.ppm[27], 401.0, atol=0.1)
    @test isapprox(res_parisext.ppm[1], 416.2, atol=0.1)
    @test isapprox(res_parisext.ppm[27], 763.5, atol=0.1)
end

@testset "Temperature increase with base and R5 scenarios" begin
    @test isapprox(res_base.TATM[1], 1.25, atol=0.051)
    @test isapprox(res_base.TATM[27], 4.91, atol=0.051)
    @test isapprox(res_r5.TATM[1], 1.25, atol=0.051)
    @test isapprox(res_r5.TATM[27], 3.24, atol=0.051)
end

@testset "Emission control rate with R4 and R3 scenarios" begin
    @test isapprox(res_r4.MIU[1], 0.05, atol=0.0051)
    @test isapprox(res_r4.MIU[17], 0.70, atol=0.0051)
    @test isapprox(res_r3.MIU[1], 0.05, atol=0.0051)
    @test isapprox(res_r3.MIU[17], 0.85, atol=0.0051)
end

@testset "Social cost of carbon with R2 and R1 scenarios" begin
    @test isapprox(res_r2.scc[1], 176, atol=0.51)
    @test isapprox(res_r2.scc[7], 302, atol=0.51)
    @test isapprox(res_r1.scc[1], 485, atol=0.51)
    @test isapprox(res_r1.scc[7], 695, atol=0.51)
end