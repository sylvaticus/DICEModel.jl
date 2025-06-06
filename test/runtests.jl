"""
Part of [DICEModel](https://github.com/sylvaticus/DICEModel.jl). Licence is MIT.

This file contains the test suit that is ran by GitHub actions on each pull.
"""


using Test, DICEModel
using JuMP, Ipopt
#using JLD2

println("Testing DICEModel...")

res_cbopt    = run_dice()
res_cbopt2   = run_dice(DICE2023())
res_cbopt3   = run_dice_scenario("cbopt")
res_t2c      = run_dice_scenario("t2c")
res_t15c     = run_dice_scenario("t15c")
res_altdam   = run_dice_scenario("altdam")
res_altdam2  = run_dice(DICE2023(a2base = [0.01]))
res_parisext = run_dice_scenario("parisext")
res_base     = run_dice_scenario("base")
res_r5       = run_dice_scenario("r5")
res_r4       = run_dice_scenario("r4")
res_r3       = run_dice_scenario("r3")
res_r2       = run_dice_scenario("r2")
res_r1       = run_dice_scenario("r1")

@testset "Test API equivalence" begin
    @test res_cbopt.ECO2 == res_cbopt2.ECO2 == res_cbopt3.ECO2
    @test res_altdam.TATM == res_altdam2.TATM
end

#for (k) in keys(res_cbopt)
#    isequal(res_cbopt[k],res_cbopt2[k]) || println(k)
#end


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

res_crazy = run_dice(a2base = 0.01, bounds = Dict("MIU_R"=>("==",1.0), "TATM"=>("<=",15), "Y_R" =>(">=",[fill(floatmin(Float64),10);fill(0.1,71)]), "ECO2" =>("<=",10000)),optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0))

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
    #@test isapprox(res_r2.scc[1], 176, atol=0.51)
    @test isapprox(res_r2.scc[7,1], 302, atol=0.51)
    #@test isapprox(res_r1.scc[1], 485, atol=0.51)
    @test isapprox(res_r1.scc[7,1], 695, atol=0.51)
end

# -----------------------------------------------------------------------------
# Scaling weigth test
@testset "Scaling weights test" begin
    w1  = [0.98,0.01,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]
    w2  = [10,20,70]
    w1b = DICEModel.scaleweights(w1)
    w2b = DICEModel.scaleweights(w2)
    @test w1b == w1
    @test sum(w2b)    == 1
    @test length(w2b) == 3
end

# -----------------------------------------------------------------------------
# Multiple regions test
res_cbopt_1r = run_dice(DICE2023_NREG(1))
res_cbopt_4r = run_dice(DICE2023_NREG(4))

@testset "DICE2023 with N equal regions" begin
    @test res_cbopt_4r.solved == true
    @test res_cbopt_1r.solved == true
    @test res_cbopt_1r.ECO2_R[:,1]  ≈ res_cbopt.ECO2_R[:,1] 
    @test res_cbopt_1r.TATM  ≈ res_cbopt.TATM 
    @test res_cbopt_4r.ECO2_R[:,1] ≈ res_cbopt_4r.ECO2_R[:,2] ≈ res_cbopt_4r.ECO2_R[:,3] ≈ res_cbopt_4r.ECO2_R[:,4]
end

w_rich  = [5,4,3,3,1,1,3,2,2,1,1.5,1]
w_equal = fill(1,12)
w_poor  = [1,1,1,1.5,2,2,3,3,2,5,5,5]

# res_cbopt_12r = run_dice(RICE2023(;weights=w_poor);optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5, "max_iter" => 1000, "acceptable_tol" =>10^-4, "acceptable_iter" => 15, "acceptable_dual_inf_tol" =>10.0^8, "acceptable_constr_viol_tol" => 0.1, "acceptable_compl_inf_tol" =>0.1, "acceptable_obj_change_tol" =>10.0^10))

res_cbopt_12r_poor = run_dice(RICE2023(;weights=w_poor)) # 657 iterations

#res_cbopt_12r_equal = run_dice(RICE2023(;weights=w_equal);optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5)) # 2006 iterations

res_cbopt_12r_rich = run_dice(RICE2023(;weights=w_rich)) # 483 iterations

@testset "RICE2023" begin
    @test res_cbopt_12r_poor.solved
    @test res_cbopt_12r_rich.solved
    # When weigth is given to the rich, US can pollute more:
    @test sum(res_cbopt_12r_rich.ECO2_R[:,1]) > sum(res_cbopt_12r_poor.ECO2_R[:,1])
end
