"""
Part of [DICEModel](https://github.com/sylvaticus/DICEModel.jl). Licence is MIT.

In this file we define the struct hosting the parameters of the model and different constructors that define each a "default" of these parameters (e.g. DICE2023 or RICE2023). User can then fine-tune this default with individual parameter choices.
"""




# -----------------------------------------------------------------------------
# Main Parameter struct (with the default of the defaults)

"""
    DICEParameters

Row and computed parameters for the optimization function.

This structure contains the "default of the default" parameters, which can eventually be modified using either keyword arguments in specific functions (each defining its own "defaults") or directly in the `run_dice(pars)` function (e.g. `run_dice(a2base = [0.01])`). This second method overrides the specific defaults of the `DICE2013` function.

The structure first defines some "raw" parameters, and then some "computed" parameters (mostly matrices of ntsteps x nregions length).
Both can be overridden with keyword arguments. In particular, "computed" parameters can be overridden in two ways: either by overriding the raw parameters from which they are computed, or by computing the parameter in a different way (outside the model) and overriding the computed parameter.

# Available parameters:
$(FIELDS)

"""
Base.@kwdef struct DICEParameters

    ######################################################################
    # Raw parameters
    ######################################################################
    
    "Years per period"
    tstep   = 5
    "Number of time periods"
    ntsteps = 81 

    "Name of the regions"
    regions = ["World"]

    "Utility weights to assign to each region. Note these are **exogenous** (default to equal weights), are NOT the Negishi weights."
    weights = [1.0/length(regions) for _ in 1:length(regions)]

    # --------------------------------------------------------------------
    # Population and technology

    "Capital elasticity in production function"
    gamma    = [0.3]    
    "Initial world population 2020 (millions)"  
    pop1    = [7752.9]
    "Growth rate to calibrate to 2050 population projection"
    popadj  =   [0.145]
    "Asymptotic population (millions)"
    popasym =  [10825]
    "Depreciation rate on capital (per year)"
    dk      = [0.1]     
    "Initial world output 2020 (trill 2019 USD)"
    q1      = [135.7]   
    "Initial level of total factor productivity"
    al1     = [5.84]    
    "Initial growth rate for TFP per 5 years"
    ga1     = [0.066]   
    "Decline rate of TFP per 5 years"
    dela    = [0.0015]  

    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output

    "Initial growth of sigma (per year)"
    gsigma1   = [-0.015]
    "Decline rate of gsigma per period"
    delgsig   = [0.96]
    "Asymptotic sigma"
    asymgsig  = [-0.005] 
    "Industrial emissions 2020 (GtCO2 per year)"
    e1        = [37.56] 
    "Emissions control rate historical 2020"
    miu1      = [0.05] 
    "Cumulative emissions 2020 (GtC)"
    cumemiss0 = 633.5

    # --------------------------------------------------------------------
    # Climate damage parameters

    "Damage intercept"
    a1     = [0]       
    "Damage quadratic term"
    a2base = [0.003467]
    "Damage exponent"
    a3     = [2]

    # --------------------------------------------------------------------
    # Abatement cost
    "Exponent of control cost function"
    expcost2  = fill(2.6,ntsteps)   
    "Cost of backstop in 2019\$ per tCO2 (2050)"
    pback2050 = [515]   

    # --------------------------------------------------------------------
    # Limits on emissions controls

    "Emission control limit from 2070"
    limmiu2070        = [1.0] 
    "Emission control limit from 2120"
    limmiu2120        = [1.1] 
    "Emission control limit from 2220"
    limmiu2200        = [1.05]
    "Emission control limit from 2300"
    limmiu2300        = [1.0] 
    "Emission control delta limit per period"
    delmiumax         = [0.12]

    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing

    "Climate beta"
    betaclim = [0.5]  
    "Elasticity of marginal utility of consumption"
    elasmu   = [0.95] 
    "Pure rate of social time preference"
    prstp    = [0.001]
    "Capital risk premium (renamed to avoid conflict with Julia's pi)"
    pi_val   = [0.05] 
    "Initial capital stock (10^12 2019 USD)"
    k0       = [295]  
    "Annual standard deviation of consumption growth"
    siggc1   = [0.01] 

    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption

    "Scaling factor for discounting"
    srf    = [1000000]   
    "Multiplicative scaling coefficient"
    scale1 = 0.00891061
    "Additive scaling coefficient"
    scale2 = -6275.91  

    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions

    "Carbon emissions from land 2015 (GtCO2 per year)"
    eland0         = [5.9]    
    "Decline rate of land emissions (per period)"
    deland         = [0.1]    
    "Non-abatable forcings 2020"
    f_misc2020     = -0.054 
    "Non-abatable forcings 2100"
    f_misc2100     = 0.265  
    "Forcings of abatable non-CO2 GHG in 2020"
    f_ghgabate2020 = [0.518]  

    "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"
    eco2eghgb2020  = [9.96]   
    "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"
    eco2eghgb2100  = [15.5]   
    "Ratio of CO2e to industrial CO2 in 2020"
    emissrat2020   = [1.4]    
    "Ratio of CO2e to industrial CO2 in 2100"
    emissrat2100   = [1.21]   
    "Coefficient of non-CO2 abateable emissions"
    fcoef1         = 0.00955
    "Coefficient of non-CO2 abateable emissions"
    fcoef2         = 0.861  

    # --------------------------------------------------------------------
    # Parameters for the DFAIR model
    
    "Calendar year that corresponds to model year zero"
    yr0      = 2020       

    "Carbon emissions share into Reservoir 0"
    emshare0 = 0.2173   
    "Carbon emissions share into Reservoir 1"
    emshare1 = 0.224    
    "Carbon emissions share into Reservoir 2"
    emshare2 = 0.2824   
    "Carbon emissions share into Reservoir 3"
    emshare3 = 0.2763   
    "Decay time constant for Reservoir 0"
    tau0     = 1000000  
    "Decay time constant for Reservoir 1"
    tau1     = 394.4    
    "Decay time constant for Reservoir 2"
    tau2     = 36.53    
    "Decay time constant for Reservoir 3"
    tau3     = 4.304    

    "Thermal equilibration parameter for box 1"
    teq1     = 0.324    
    "Thermal equilibration parameter for box 2"
    teq2     = 0.44     
    "Thermal response timescale for deep ocean"
    d1       = 236      
    "Thermal response timescale for upper ocean"
    d2       = 4.07     

    "Pre-industrial IRF100"
    irf0     = 32.4     
    "Increase in IRF100 with cumulative carbon uptake"
    irc      = 0.019    
    "Increase in IRF100 with warming"
    irt      = 4.165    
    "Forcings of equilibrium CO2 doubling"
    fco22x   = 3.93     
    # --------------------------------------------------------------------
    # Initial conditions to be calibrated to history calibration

    "Initial concentration in atmosphere in 2020 (GtC)"
    mat0   = 886.5128014

    "Initial concentration in Reservoir 0 in 2020 (GtC)"
    res00  = 150.093    
    "Initial concentration in Reservoir 1 in 2020 (GtC)"
    res10  = 102.698    
    "Initial concentration in Reservoir 2 in 2020 (GtC)"
    res20  = 39.534     
    "Initial concentration in Reservoir 3 in 2020 (GtC)"
    res30  = 6.1865     

    "Equilibrium concentration in atmosphere (GtC)"
    mateq  = 588        
    "Initial temperature box 1 change in 2020 (°C)"
    tbox10 = 0.1477     
    "Initial temperature box 2 change in 2020 (°C)"
    tbox20 = 1.099454   
    "Initial atmospheric temperature change in 2020 (°C)"
    tatm0  = 1.247154    # Changed to prevent numerical instability (GAMS original: 1.24715)

    ######################################################################
    # Computed parameters
    ######################################################################

    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing

    "Time periods sequence (0,5,10,...,400)"
    times= 0:tstep:(ntsteps*tstep-1)
    "Time periods index sequence (1,2,3,...,81)"
    tidx  = 1:ntsteps           
    "Time periods index sequence (0,1,2,...,80)"      
    t0idx = 0:ntsteps-1    
    "Number of regions"
    nreg  = length(regions)
    "Regions index sequence"
    ridx = 1:nreg

    "Risk-adjusted rate of time preference"
    rartp = @. exp(prstp + betaclim*pi_val)-1  

    # --------------------------------------------------------------------
    # Limits on emissions controls (computed)

    "Upper bounds on miu"
    miuup = [
        if ti == 1
            0.05
        elseif ti == 2
            0.10
        elseif ti < 9
            delmiumax[ri]*(ti - 1)
        elseif ti < 12
            0.85+.05*(ti-8)
        elseif ti < 21
            limmiu2070[ri]
        elseif ti < 38
            limmiu2120[ri]
        elseif ti < 58
            limmiu2200[ri]
        else
            limmiu2300[ri]       
        end
        for ti in tidx, ri in ridx
    ]

    # --------------------------------------------------------------------
    # Precautionary parameters

    "Variance of per capita consumption"
    varpcc    = [min((siggc1[ri]^2)*t,(siggc1[ri]^2)*tstep*47) for t in times, ri in ridx]
    "Precautionary rate of return" 
    rprecaut  = [-0.5 * varpcc[ti,ri] * elasmu[ri]^2 for ti in tidx, ri in ridx]
    "STP factor without precautionary factor"
    rr1       = [1 / ((1+rartp[ri])^times[ti]) for ti in tidx, ri in ridx]
    "STP with precautionary factor"
    rr        = [rr1[ti,ri] * (1+rprecaut[ti,ri]) ^ -times[ti] for ti in tidx, ri in ridx]
    "Optimal long-run savings rate used for transversality"
    optlrsav  = @. (dk + 0.004)/(dk + 0.004*elasmu + rartp)*gamma

    # --------------------------------------------------------------------
    # Dynamic parameters

    "Level of population and labor (temp, used only for its first value)"
    l_temp     = fill(0.0, ntsteps,nreg)
    "Level of population and labor"
    l          = [l_temp[ti,ri] = (ti == 1) ? pop1[ri] : l_temp[ti-1,ri] * (popasym[ri] / l_temp[ti-1,ri])^popadj[ri] for ti in tidx, ri in ridx]
    "Growth rate of Total Factor Productivity"
    ga         = [ga1[ri]*exp(-dela[ri]*times[ti]) for ti in tidx, ri in ridx]
    "Level of total factor productivity (temp, used only for its first value)"
    al_temp    = fill(0.0,ntsteps,nreg)
    "Level of total factor productivity"
    al         = [al_temp[ti,ri] = (ti == 1) ? al1[ri] :  al_temp[ti-1,ri] / (1 - ga[ti-1,ri]) for ti in tidx,ri in ridx]

    "Backstop price 2019\$ per ton CO2"
    pbacktime  = hcat([vcat(pback2050[ri] .* exp.(-tstep .* 0.01 .* (tidx[1:7] .- 7)), pback2050[ri] .* exp.(-tstep .* 0.001 .*(tidx[8:end] .-7 )) ) for ri in ridx]...)

    "Carbon intensity 2020 kgCO2-output 2020"
    sig1       = @. e1/(q1*(1-miu1))

    "Change in sigma (rate of decarbonization)"
    gsig       = [min(gsigma1[ri]*delgsig[ri] ^(ti-1),asymgsig[ri]) for ti in tidx, ri in ridx]

    "CO2-emissions output ratio (temp, used only for its first value)"
    sigma_temp = fill(0.0,ntsteps, nreg)
    "CO2-emissions output ratio"
    sigma      = [sigma_temp[ti,ri] = (ti==1) ? sig1[ri] : sigma_temp[ti-1,ri] * exp(tstep*gsig[ti-1,ri]) for ti in tidx, ri in ridx]

    # ------------------------------------------------------------------------------
    # Parameters emissions and non-CO2 
    
    eland          = [eland0[ri]*(1-deland[ri])^ti for ti in t0idx, ri in ridx]      
    co2e_ghgabateb = [eco2eghgb2020[ri] + (eco2eghgb2100[ri]-eco2eghgb2020[ri]) * min(1,ti/16) for ti in t0idx, ri in ridx]
    f_misc         = f_misc2020    .+ [(f_misc2100-f_misc2020) * min(1,ti/16) for ti in t0idx]
    emissrat       = [emissrat2020[ri] + (emissrat2100[ri]-emissrat2020[ri]) * min(1,ti/16) for ti in t0idx, ri in ridx]
    sigmatot       = @. sigma * emissrat
    cost1tot       = [pbacktime[ti,ri] * sigmatot[ti,ri] / expcost2[ti,ri] / 1000 for ti in tidx, ri in ridx]

end

# -----------------------------------------------------------------------------
# Different constructors for the different defaults

"""
  DICE2023()

Parameters constructor with defaults aligned to DICE2023. 

Create a parameters struct with the defaults of the DICE2023 model. Different parameters, either the "raw" ones or the "computed ones", can be specified using keywork arguments.

See [`DICEParameters`](@ref) for a complete list of available parameters.

# Example
```julia
alt_dam_scen = DICE2023(a2base = [0.01]) # single value array parameter for the "world region"
```

# Note
- Even if DICE2023 treats the World as a single region, the model works with a "regional" dimension, so all parameters (except those related to the carbon cycle) should be entered as a single-value array for static data or a (ntime_periods by nregions) matrix for (computed) dynamic parameters.

"""
function DICE2023(;kwargs...)
    pars = DICEParameters(;kwargs...) 
    return pars
end

"""
    DICE2023_NREG(n;kwargs...)

Build parameters for a DICE2023 world partitioned in n equal regions (unless parameters are overrided)

Create a parameters struct where the world is partitioned in n regions, and where each region is equal, with the same coefficients but with 1/n of initial emissions, capital, production, population, ...

Data is from DICE2013, so the output of `DICE2023_NREG(1)` is the same as `DICE2023()`.

Using the keyword arguments one can specify individual parameters that differ from DICE2023, with eventually a regional specification.
    
See [`DICEParameters`](@ref) for a complete list of available parameters and [`run_dice`](@ref) to run the optimization with the parameter struct created with this function.

# Example
```julia
four_regions_parameters = DICE2023_NREG(4) 
alt_dam_parameters      = DICE2023_NREG(2, a2base = [0.01, 0.02]) # two regions differing only for the a2base parameter
results = run_dice(four_regions_parameters)
```
"""
function DICE2023_NREG(n=2;
    regions = ["Reg_$(r_n)" for r_n in 1:n],
    tstep   = 5, # "Years per period"
    ntsteps = 81, # "Number of time periods"
    weights = fill(1/n,n),
    gamma   = fill(0.3,n), # "Capital elasticity in production function"    
    pop1    = fill(7752.9/n,n), # "Initial world population 2020 (millions)" 
    popadj  = fill(0.145,n), # "Growth rate to calibrate to 2050 population projection"
    popasym = fill(10825/n,n), # "Asymptotic population (millions)"
    dk      = fill(0.1,n), # "Depreciation rate on capital (per year)"     
    q1      = fill(135.7/n,n), # "Initial world output 2020 (trill 2019 USD)" 
    al1     = fill(5.84,n), # "Initial level of total factor productivity"    
    ga1     = fill(0.066,n), # "Initial growth rate for TFP per 5 years"   
    dela    = fill(0.0015,n), # "Decline rate of TFP per 5 years"
    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    gsigma1   = fill(-0.015,n), # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,n), # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,n), # "Asymptotic sigma"
    e1        = fill(37.56/n,n), # "Industrial emissions 2020 (GtCO2 per year)" 
    miu1      = fill(0.05,n), # "Emissions control rate historical 2020" 
    # --------------------------------------------------------------------
    # Climate damage parameters    
    a1     = fill(0,n), # "Damage intercept"       
    a2base = fill(0.003467,n), # "Damage quadratic term"
    a3     = fill(2,n), # "Damage exponent"
    # --------------------------------------------------------------------
    # Abatement cost
    expcost2  = fill(2.6,ntsteps,n), # "Exponent of control cost function"
    pback2050 = fill(515,n), # "Cost of backstop in 2019\$ per tCO2 (2050)"
    # --------------------------------------------------------------------
    # Limits on emissions controls
    limmiu2070        = fill(1.0,n), # "Emission control limit from 2070"
    limmiu2120        = fill(1.1,n), # "Emission control limit from 2120"
    limmiu2200        = fill(1.05,n), # "Emission control limit from 2220"
    limmiu2300        = fill(1.0,n), # "Emission control limit from 2300"
    delmiumax         = fill(0.12,n), # "Emission control delta limit per period"
    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing
    betaclim = fill(0.5,n), # "Climate beta"  
    elasmu   = fill(0.95,n), # "Elasticity of marginal utility of consumption" 
    prstp    = fill(0.001,n), # "Pure rate of social time preference"
    pi_val   = fill(0.05,n), # "Capital risk premium (renamed to avoid conflict with Julia's pi)" 
    k0       = fill(295/n,n), # "Initial capital stock (10^12 2019 USD)"  
    siggc1   = fill(0.01,n), # "Annual standard deviation of consumption growth"
    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    srf    = fill(1000000,n), # "Scaling factor for discounting"
    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions  
    eland0         = fill(5.9/n,n), # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = fill(0.1,n), # "Decline rate of land emissions (per period)"    
    f_ghgabate2020 = fill(0.518/n,n), # "Forcings of abatable non-CO2 GHG in 2020"  
    eco2eghgb2020  = fill(9.96/n,n), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"     
    eco2eghgb2100  = fill(15.5/n,n), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"   
    emissrat2020   = fill(1.4,n), # "Ratio of CO2e to industrial CO2 in 2020"    
    emissrat2100   = fill(1.21,n), # "Ratio of CO2e to industrial CO2 in 2100"
    kwargs...
    )
    pars = DICEParameters(;regions, tstep, ntsteps, weights, gamma, pop1, popadj, popasym, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, eland0, deland, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, kwargs...) 
    return pars
end

"""
    RICE2023(;kwargs...)

Build parameters calibrated to have the DICE2023 world totals and the 12-regions RICE2020 regional distribution.

Create a parameters struct where the world is partitioned in 12 regions, with coefficients of DICE2023 but a regional variance derived in most cases from RICE2020 (initial emissions, capital, production, population) or assumed (e.g. damages)

Note that the utility weights to provide to each region are exogenous (default to equal weights), they are NOT the Negishi weights. The `run_dice` function can eventually be used iteractively to look for these weights.

See [`DICEParameters`](@ref) for a complete list of available parameters and [`run_dice`](@ref) to run the optimization with the parameter struct created with this function.

# Example
```julia
w_rich  = [5,4,3,3,1,1,3,2,2,1,1.5,1] #
w_equal = fill(1,12)
res_cbopt_12r_poor = run_dice(RICE2023(;weights=w_poor))
res_cbopt_12r_rich = run_dice(RICE2023(;weights=w_rich)) 
```
# Notes
- The default 12 regions are: ["USA", "EUS", "JPN", "OHI", "RUS", "EEC", "CHN", "IND", "MDE", "SSA", "LAA", "ROW"]
"""
function RICE2023(;
    regions = ["USA", "EUS", "JPN", "OHI", "RUS", "EEC", "CHN", "IND", "MDE", "SSA", "LAA", "ROW"],
    weights = fill(1,12),
    tstep   = 5, # "Years per period"
    ntsteps = 81, # "Number of time periods"
    # ----------------------------------
    # Main regionalised data
    
    # Population
    pop1    = [321, 320, 128, 244, 145, 149, 1407, 1310, 407, 930, 615, 1369], # "Initial world population 2020 (millions)" 
    popadj  = [0.08, 0.03, 0.03, 0.06, 0.03, 0.03, 0.03, 0.2, 0.15, 0.2, 0.2, 0.2],     # "Growth rate to calibrate to 2050 population projection"
    popasym = [500, 360, 160, 370, 170, 170, 1500, 2100, 600, 1800, 1000, 2100],  # "Asymptotic population (millions)"
    # Production and factor productivity
    q1      = [21.7, 16.4, 6.1, 12.0, 4.6, 3.5, 23.6, 9.6, 10.4, 4.2, 10.7, 13.1],  # "Initial world output 2020 (trill 2019 USD)" 
    k0      = [43.3, 37.8, 16.8, 27.9, 7.3, 6.3, 64.7, 18.0, 18.2, 7.3, 22.8, 24.6],    # "Initial capital stock (10^12 2019 USD)"  
    al1     = fill(5.84,12),      # "Initial level of total factor productivity"    
    ga1     = fill(0.066,12),     # "Initial growth rate for TFP per 5 years"   
    dela    = fill(0.0015,12),    # "Decline rate of TFP per 5 years"
    gamma   = fill(0.3,12),       # "Capital elasticity in production function"
    # Emission intensity
    e1        = [5.71, 2.40, 1.30, 2.62, 1.89, 0.93, 11.30, 2.61, 3.01, 0.89, 1.98, 2.92], # "Industrial emissions 2020 (GtCO2 per year)" 
    gsigma1   = fill(-0.015,12),   # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,12),     # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,12),   # "Asymptotic sigma"
    # Emission from land use
    eland0         = [0,0,0,0,0,0,0.4,0,0,0,0,5.5], # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = fill(0.1,12),    # "Decline rate of land emissions (per period)"  
    # Cost of mitigation
    # The regional diversity is already accounted in the computation of cost1tot variable, based on emission intensity of the different economies
    expcost2  = fill(2.6,ntsteps,12),     # "Exponent of control cost function"

    # Climate change damage
    # 2 region example total dam equivalence with dam in 2nd region double (in proportion to Y terms): 
    # a2base*T^a3*Y == a2base_1 * T^a3*Y_1 + 2 * a2base_1 * T^a3 * Y_2
    # ==> a2base_1 = a2base * Y/(Y_1 + 2 Y_2)
    # More in general: a2_i = a2 * sum_j(Y_j) / sum_j (coef_j/coef_i Y_i))
    # dam coefficients (d):  [1, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 2, 2, 2]
    # A2 base: 0.003467
    # a2base_r = [0.003467 * sum(q1)/ sum((d[j]/d[i]) * q1[j] for j in ridx ) for i in ridx]
    a2base = [0.002709, 0.002709, 0.002709, 0.002709, 0.002709, 0.002709, 0.002709, 0.004064, 0.004064, 0.005419, 0.005419, 0.005419], # "Damage quadratic term"
    # --------------------------------------------------------------------     
    dk      = fill(0.1,12), # "Depreciation rate on capital (per year)"   
    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    
    miu1      = fill(0.05,12), # "Emissions control rate historical 2020" 
    # --------------------------------------------------------------------
    # Climate damage parameters    
    a1     = fill(0,12), # "Damage intercept"       
    a3     = fill(2,12), # "Damage exponent"
    # --------------------------------------------------------------------
    # Abatement cost
    pback2050 = fill(515,12), # "Cost of backstop in 2019\$ per tCO2 (2050)"
    # --------------------------------------------------------------------
    # Limits on emissions controls
    limmiu2070        = fill(1.0,12), # "Emission control limit from 2070"
    limmiu2120        = fill(1.1,12), # "Emission control limit from 2120"
    limmiu2200        = fill(1.05,12), # "Emission control limit from 2220"
    limmiu2300        = fill(1.0,12), # "Emission control limit from 2300"
    delmiumax         = fill(0.12,12), # "Emission control delta limit per period"
    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing
    betaclim = fill(0.5,12), # "Climate beta"  
    elasmu   = fill(0.95,12), # "Elasticity of marginal utility of consumption" 
    prstp    = fill(0.001,12), # "Pure rate of social time preference"
    pi_val   = fill(0.05,12), # "Capital risk premium (renamed to avoid conflict with Julia's pi)" 
    siggc1   = fill(0.01,12), # "Annual standard deviation of consumption growth"
    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    srf    = fill(1000000,12), # "Scaling factor for discounting"
    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions  
    f_ghgabate2020 = fill(0.518/12,12), # "Forcings of abatable non-CO2 GHG in 2020"  
    eco2eghgb2020  = fill(9.96/12,12), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"     
    eco2eghgb2100  = fill(15.5/12,12), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"   
    emissrat2020   = fill(1.4,12), # "Ratio of CO2e to industrial CO2 in 2020"    
    emissrat2100   = fill(1.21,12), # "Ratio of CO2e to industrial CO2 in 2100"
    kwargs...
    )
    pars = DICEParameters(;regions, tstep, ntsteps,weights, gamma, pop1, popadj, popasym, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, eland0, deland, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, kwargs...) 
    return pars
end
