"""
    run_dice_scenario(scenario::String)

Run one of the "official" 11 scenarios in Nordhous's DICE 2023 model:

- `cbopt`: The C/B optimal scenario
- `t2c`:   The temperature constrained to max 2 °C scenario
- `t15c`:   The temperature constrained to max 1.5 °C scenario
- `altdam`: The alternative damage scenario
- `parisext`: The Paris extended scenario
- `base`:     The base (current policies) scenario
- `r5`:       The scenario with discount rate of 5%
- `r4`:       The scenario with discount rate of 4%
- `r3`:       The scenario with discount rate of 3%
- `r2`:       The scenario with discount rate of 2%
- `r1`:       The scenario with discount rate of 1%

To run "your own" scenarios, use the function [`run_dice`](@ref), where you can set the input parameters and constraints as you like.

# Output

- A named tuple containing the following fields: `solved`, `status`, `times`, `tidx`, the post_process computed values and the optimization variables.

"""
function run_dice_scenario(scenario::String)

    if !(scenario in ["cbopt","t2c","t15c","altdam","parisext","base","r5","r4","r3","r2","r1"])
        error("Unknown scenario $scenario")
    end

    if scenario == "cbopt"
        return run_dice() # run_dice is already by default the optimal CB scenario
    elseif scenario == "t2c"
        return run_dice(bounds = Dict("TATM"=>("<=",2.0)))
    elseif scenario == "t15c"
        return run_dice(bounds = Dict("TATM"=>("<=",1.5)))
    elseif scenario == "altdam"
        return run_dice(a2base = 0.01)
    elseif scenario == "parisext"
        tidx = 1:81
        miuup = @. min( 0.05 + 0.04*(tidx-1) - 0.01*max(0,tidx-5)  ,1.00)
        return run_dice(miuup = miuup)
    elseif scenario == "base"
        times   = 0:5:400
        cprice1 = 6      # Carbon price in 2020 (2019$ per tCO2)
        gcprice = 0.025  # Growth rate of base carbon price per year
        miuup   = fill(1.0,81)
        
        # Carbon price in base case
        cpricebase  = @. cprice1*(1+gcprice)^times
        cpriceub    = [cpricebase[1:46]; fill(floatmax(Float64),81-46)]
        miulb = [fill(0.0,57);fill(1.0,81-57)]
        bounds = Dict("CPRICE_R"=>("<=",cpriceub),"MIU_R"=>(">=",miulb))
        return run_dice(miuup = miuup, bounds = bounds)
    elseif scenario == "r5"
        prstp     = 0.05
        elasmu    = 0.001
        tidx      = 1:81
        times     = 0:5:400
        dk        = 0.1 
        gama      = 0.3 
        rr1       = @. 1/(1+prstp)^times
        [rr1[tid] = 1/(1+prstp)^(5*51) for tid in tidx[52:end]]
        rr        = rr1
        optlrsav  = (dk + 0.004)/(dk + 0.004*elasmu + prstp)*gama
        k0        = 290
        return run_dice(prstp = prstp, elasmu = elasmu,rr1=rr1,rr=rr,k0=k0,optlrsav = optlrsav)
    elseif scenario == "r4"
        prstp     = 0.04
        elasmu    = 0.001
        tidx      = 1:81
        times     = 0:5:400
        dk        = 0.1 
        gama      = 0.3 
        rr1       = @. 1/(1+prstp)^times
        [rr1[tid] = 1/(1+prstp)^(5*80) for tid in tidx[82:end]]
        rr        = rr1
        optlrsav  = (dk + 0.004)/(dk + 0.004*elasmu + prstp)*gama
        k0        = 326
        return run_dice(prstp = prstp, elasmu = elasmu,rr1=rr1,rr=rr,k0=k0,optlrsav = optlrsav,bounds=Dict("CPRICE_R"=>("<=",1000)))
    elseif scenario == "r3"
        prstp     = 0.03
        elasmu    = 0.001
        tidx      = 1:81
        times     = 0:5:400
        dk        = 0.1 
        gama      = 0.3 
        rr1       = @. 1/(1+prstp)^times
        rr        = rr1
        optlrsav  = (dk + 0.004)/(dk + 0.004*elasmu + prstp)*gama
        k0        = 370
        return run_dice(prstp = prstp, elasmu = elasmu,rr1=rr1,rr=rr,k0=k0,optlrsav = optlrsav)
    elseif scenario == "r2"
        prstp     = 0.02
        elasmu    = 0.001
        tidx      = 1:81
        times     = 0:5:400
        dk        = 0.1 
        gama      = 0.3 
        rr1       = @. 1/(1+prstp)^times
        rr        = rr1
        optlrsav  = (dk + 0.004)/(dk + 0.004*elasmu + prstp)*gama
        k0        = 409
        return run_dice(prstp = prstp, elasmu = elasmu,rr1=rr1,rr=rr,k0=k0,optlrsav = optlrsav)
    elseif scenario == "r1"
        prstp     = 0.01
        elasmu    = 0.001
        tidx      = 1:81
        times     = 0:5:400
        dk        = 0.1 
        gama      = 0.3 
        rr1       = @. 1/(1+prstp)^times
        rr        = rr1
        optlrsav  = (dk + 0.004)/(dk + 0.004*elasmu + prstp)*gama
        k0        = 420
        return run_dice(prstp = prstp, elasmu = elasmu,rr1=rr1,rr=rr,k0=k0,optlrsav = optlrsav)
    end
end