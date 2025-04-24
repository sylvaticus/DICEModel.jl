"""
In this file we define the struct hosting the parameters of the model and different constructors that define each a "default" of these parameters (e.g. DICE2023 or RICE2023). User can then fine-tune this default with individual parameter choices.
"""




# -----------------------------------------------------------------------------
# Main Parameter struct (with the default of the defaults)

"""
    DICEParameters

This structure contains the "default" parameters, which can eventually be modified using keyword arguments in the `run_dice(pars)` function (e.g. `run_dice(a2base = [0.01])`).

The structure first defines some "raw" parameters, and then some "computed" parameters (mostly arrays of ntsteps length).
Both can be overridden with keyword arguments in the `run_dice(pars)` function. In particular, "computed" parameters can be overridden in two ways: either by overriding the raw parameters from which they are computed, or by computing the parameter in a different way (outside the model) and overriding the computed parameter.

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

    "Weights"
    weights = [1.0]

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
    expcost2  = [2.6]   
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
    cost1tot       = [pbacktime[ti,ri] * sigmatot[ti,ri] / expcost2[ri] / 1000 for ti in tidx, ri in ridx]

end

# -----------------------------------------------------------------------------
# Different constructors for the different defaults

"""
  DICE2023()

Parameters constructor with defaults aligned to DICE2023. 

Create a parameters struct with the defaults of DICE2023 model. Different parameters, either the "raw" ones or the "computed ones", can be specified using keywork arguments.

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
  RICE2020()

Parameters constructor with defaults aligned to RICE2020. 

Create a parameters struct with the defaults of the 12-regions RICE2020 model. Different parameters, either the "raw" ones or the "computed ones", can be specified using keywork arguments.

See [`DICEParameters`](@ref) for a complete list of available parameters.

# Example
```julia
alt_dam_scen = DICEParameters(a2base = [0.01]) # single value array parameter for the "world region"
```

# Note
- Something in the middle struct of DICE2023 and data of RICE2020

"""
function DICE2023_2REG(;
    regions = ["World_half1","World_half2"],

    weights = fill(0.5,2),
    gamma    = fill(0.3,2), # "Capital elasticity in production function"    
    pop1    = fill(7752.9/2,2), # "Initial world population 2020 (millions)" 
    popadj  = fill(0.145,2), # "Growth rate to calibrate to 2050 population projection"
    popasym = fill(10825/2,2), # "Asymptotic population (millions)"
    dk      = fill(0.1,2), # "Depreciation rate on capital (per year)"     
    q1      = fill(135.7/2,2), # "Initial world output 2020 (trill 2019 USD)" 
    al1     = fill(5.84,2), # "Initial level of total factor productivity"    
    ga1     = fill(0.066,2), # "Initial growth rate for TFP per 5 years"   
    dela    = fill(0.0015,2), # "Decline rate of TFP per 5 years"
    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    gsigma1   = fill(-0.015,2), # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,2), # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,2), # "Asymptotic sigma"
    e1        = fill(37.56/2,2), # "Industrial emissions 2020 (GtCO2 per year)" 
    miu1      = fill(0.05,2), # "Emissions control rate historical 2020" 
    # --------------------------------------------------------------------
    # Climate damage parameters    
    a1     = fill(0,2), # "Damage intercept"       
    a2base = fill(0.003467,2), # "Damage quadratic term"
    a3     = fill(2,2), # "Damage exponent"
    # --------------------------------------------------------------------
    # Abatement cost
    expcost2  = fill(2.6,2), # "Exponent of control cost function"
    pback2050 = fill(515,2), # "Cost of backstop in 2019\$ per tCO2 (2050)"
    # --------------------------------------------------------------------
    # Limits on emissions controls
    limmiu2070        = fill(1.0,2), # "Emission control limit from 2070"
    limmiu2120        = fill(1.1,2), # "Emission control limit from 2120"
    limmiu2200        = fill(1.05,2), # "Emission control limit from 2220"
    limmiu2300        = fill(1.0,2), # "Emission control limit from 2300"
    delmiumax         = fill(0.12,2), # "Emission control delta limit per period"
    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing
    betaclim = fill(0.5,2), # "Climate beta"  
    elasmu   = fill(0.95,2), # "Elasticity of marginal utility of consumption" 
    prstp    = fill(0.001,2), # "Pure rate of social time preference"
    pi_val   = fill(0.05,2), # "Capital risk premium (renamed to avoid conflict with Julia's pi)" 
    k0       = fill(295/2,2), # "Initial capital stock (10^12 2019 USD)"  
    siggc1   = fill(0.01,2), # "Annual standard deviation of consumption growth"
    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    srf    = fill(1000000,2), # "Scaling factor for discounting"
    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions  
    eland0         = fill(5.9/2,2), # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = fill(0.1,2), # "Decline rate of land emissions (per period)"    
    f_ghgabate2020 = fill(0.518/2,2), # "Forcings of abatable non-CO2 GHG in 2020"  
    eco2eghgb2020  = fill(9.96/2,2), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"     
    eco2eghgb2100  = fill(15.5/2,2), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"   
    emissrat2020   = fill(1.4,2), # "Ratio of CO2e to industrial CO2 in 2020"    
    emissrat2100   = fill(1.21,2), # "Ratio of CO2e to industrial CO2 in 2100"
    kwargs...
    )
    pars = DICEParameters(;regions, weights, gamma, pop1, popadj, popasym, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, eland0, deland, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, kwargs...) 
    return pars
end

function RICE2020(;
    tstep   = 5, # "Years per period"
    ntsteps = 81, # "Number of time periods" 
    regions = ["USA", "EUS", "JPN", "OHI", "RUS", "EEC", "CHN", "IND", "MDE", "SSA", "LAA", "ROW"],
    y0      = [17.060, 12.898, 4.810, 9.452, 3.609, 2.720, 18.559, 7.525, 8.157, 3.338, 8.434, 10.348],  # Y0 "Initial world output 2015 (trill 2019 USD)"    
    k0      = [32.957, 28.793, 12.777, 21.242, 5.593, 4.830, 49.296, 13.684, 13.820, 5.562, 17.346, 18.717], # K0 "Initial capital stock (10^12 2019 USD)" 
    l0      = [0.321, 0.320, 0.128, 0.244, 0.145, 0.149, 1.407, 1.310, 0.407, 0.930, 0.615, 1.369] .* 1000,
    e_f0    = [1.398, 0.588, 0.317, 0.642, 0.463, 0.228, 2.767, 0.638, 0.736, 0.219, 0.486, 0.715],
    e_l0    = [0, 0, 0, 0, 0, 0, 0.04, 0, 0, 0, 0, 0.56],

    lgr     = [3.73 , 1.08 , 4.28 , 4.28 , 2.08 , 2.08 , 3.07 , 8.09 , 8.09 , 8.09 , 8.09 , 8.09],
    lgrgr   = [11.12 , 101.33 , 50.88 , 50.88 , 67.36 , 67.36 , 97.13 , 12.44 , 12.44 , 12.44 , 12.44 , 12.44],
    tfpgr   = [0.045, 0.065, 0.065, 0.065, 0.08, 0.08, 0.08, 0.07, 0.07, 0.07, 0.07, 0.07],
    tfpgrgr = [0.0095, 0.011, 0.011, 0.011, 0.01, 0.01, 0.00965, 0.009, 0.009, 0.009, 0.009, 0.009],
    sig_i   = [9, 10, 12, 12, 10, 7.5, 8, 8, 7, 7, 7, 7],
    sig_a   = [0.14, 0.2, 0.25, 0.25, 0.2, 0.15, 0.2, 0.2, 0.15, 0.15, 0.15, 0.15],

    a1_fp   = [0.01102, 0.02, 0.01174, 0.01174, 0.02093, 0.01523, 0.02, 0.02093, 0.021, 0.021, 0.02093, 0.02093], # had to rename as a1 already used in DICE2023 for other stuff
    a2_fp   = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5], # renamed for consistency with a1_fp
    b1      = [0.07, 0.1, 0.05, 0.05, 0.1, 0.15, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
    b2      = [2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887, 2.887],
    dam1    = [-0.0026, -0.001, -0.0042, -0.0042, -0.001, 0.0108, -0.0041, -0.0041, 0.0039, 0.0039, 0.0039, 0.0039],
    dam2    = [0.0017, 0.0049, 0.0025, 0.0025, 0.0049, 0.0033, 0.002, 0.002, 0.0013, 0.0013, 0.0013, 0.0013],

    nreg    = length(regions),
    weights = fill(1/nreg,nreg),
    gamma   = fill(0.3,nreg), # "Capital elasticity in production function"     

    dk      = fill(0.1,nreg), # "Depreciation rate on capital (per year)"       
    al1     = fill(5.84,nreg), # "Initial level of total factor productivity"    
    ga1     = fill(0.066,nreg), # "Initial growth rate for TFP per 5 years"   
    dela    = fill(0.0015,nreg), # "Decline rate of TFP per 5 years"
    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    gsigma1   = fill(-0.015,nreg), # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,nreg), # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,nreg), # "Asymptotic sigma"
    e1        = fill(37.56/2,nreg), # "Industrial emissions 2020 (GtCO2 per year)" 
    miu1      = fill(0.05,nreg), # "Emissions control rate historical 2020" 
    # --------------------------------------------------------------------
    # Climate damage parameters    
    a1     = fill(0,nreg), # "Damage intercept"       
    a2base = fill(0.003467,nreg), # "Damage quadratic term"
    a3     = fill(2,nreg), # "Damage exponent"
    # --------------------------------------------------------------------
    # Abatement cost
    expcost2  = fill(2.6,nreg), # "Exponent of control cost function"
    pback2050 = fill(515,nreg), # "Cost of backstop in 2019\$ per tCO2 (2050)"
    # --------------------------------------------------------------------
    # Limits on emissions controls
    limmiu2070        = fill(1.0,nreg), # "Emission control limit from 2070"
    limmiu2120        = fill(1.1,nreg), # "Emission control limit from 2120"
    limmiu2200        = fill(1.05,nreg), # "Emission control limit from 2220"
    limmiu2300        = fill(1.0,nreg), # "Emission control limit from 2300"
    delmiumax         = fill(0.12,nreg), # "Emission control delta limit per period"

    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing
    betaclim = fill(0.5,nreg), # "Climate beta"  
    elasmu   = fill(0.95,nreg), # "Elasticity of marginal utility of consumption" 
    prstp    = fill(0.001,nreg), # "Pure rate of social time preference"
    pi_val   = fill(0.05,nreg), # "Capital risk premium (renamed to avoid conflict with Julia's pi)" 
    siggc1   = fill(0.01,nreg), # "Annual standard deviation of consumption growth"
    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    srf    = fill(1000000,nreg), # "Scaling factor for discounting"
    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions    
    f_ghgabate2020 = fill(0.518/nreg,nreg), # "Forcings of abatable non-CO2 GHG in 2020"  
    eco2eghgb2020  = fill(9.96/nreg,nreg), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"     
    eco2eghgb2100  = fill(15.5/nreg,nreg), # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"   
    emissrat2020   = fill(1.4,nreg), # "Ratio of CO2e to industrial CO2 in 2020"    
    emissrat2100   = fill(1.21,nreg), # "Ratio of CO2e to industrial CO2 in 2100"
    kwargs...
    )
   
    times = 0:tstep:(ntsteps*tstep-1) #  "Time periods sequence (0,5,10,...,400)"
    tidx  = 1:ntsteps           # "Time periods index sequence (1,2,3,...,81)"
    t0idx = 0:ntsteps-1    # "Time periods index sequence (0,1,2,...,80)"      
    ridx = 1:nreg # "Regions index sequence"
    
    # Parameter computations that differ from DICE 2023...
    # Population
    pop1 = l0
    lcgr = [lgr[r]/lgrgr[r] * (1-exp(- (t-1)* lgrgr[r]/100)) for t in tidx, r in ridx] 
    l    = [pop1[r] * exp(lcgr[t,r]) for t in tidx, r in ridx]

    # Production and productivity
    q1 = y0 
    tfp0 =  [q1[r]*(1+a1_fp[r]*(0.8/2.5)^a2_fp[r]) / ((k0[r]^gamma[r]) * (pop1[r]^(1-gamma[r]))) for r in ridx]
    al  =  [tfp0[r]*exp(tfpgr[r] * exp(-tfpgrgr[r]*t) * t) for t in t0idx, r in ridx] # total factor productivity

    # Carbon intensity
    ssig  = log.(sig_a)
    sig0  = e_f0 ./ q1
    dsig  = @. log(1-log(1+sig_i/100) / ssig )
    gsig  = [ssig[r] * (1-exp(-dsig[r]*t)) for t in t0idx, r in ridx]
    sigma = [sig0[r] * exp(gsig[t,r]) for t in tidx, r in ridx] # carbon intensity

    # Non-industrial non-abatable emissions
    eland = [e_l0[r] * (1-0.1)^(t/2) for t in tidx, r in ridx]

    # Damages
    # Variable: a2base

    # Cost of mitigation
    # Variable: expcost2


    pars = DICEParameters(;regions, weights, gamma, pop1, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, l, al, sigma, eland, kwargs...) 
    return pars
end

"""
  RICE2023

Most regionalised data is the data in RICE2020 rescaled to match DICE2023 totals.

"""
function RICE2023(;
    regions = ["USA", "EUS", "JPN", "OHI", "RUS", "EEC", "CHN", "IND", "MDE", "SSA", "LAA", "ROW"],
    weights = fill(1,12),
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
    gsigma1   = fill(-0.015,12),   # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,12),     # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,12),   # "Asymptotic sigma"
    e1        = [5.71, 2.40, 1.30, 2.62, 1.89, 0.93, 11.30, 2.61, 3.01, 0.89, 1.98, 2.92], # "Industrial emissions 2020 (GtCO2 per year)" 
    # Emission from land use
    eland0         = [0,0,0,0,0,0,0.4,0,0,0,0,5.5], # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = fill(0.1,12),    # "Decline rate of land emissions (per period)"  

    # Cost of mitigation
    expcost2  = fill(2.6,12),     # "Exponent of control cost function"

    # Climate change damage
    a2base = fill(0.003467,12), # "Damage quadratic term"

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
    pars = DICEParameters(;regions, weights, gamma, pop1, popadj, popasym, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, eland0, deland, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, kwargs...) 
    return pars
end



#=
function RICE2023(;
    regions = ["USA", "EUS", "JPN", "OHI", "RUS", "EEC", "CHN", "IND", "MDE", "SSA", "LAA", "ROW"],
    weights = fill(1,12),
    # ----------------------------------
    # Main regionalised data
    
    # Population
    pop1    = fill(7752.9/12,12), # "Initial world population 2020 (millions)" 
    popadj  = fill(0.145,12),     # "Growth rate to calibrate to 2050 population projection"
    popasym = fill(10825/12,12),  # "Asymptotic population (millions)"
    # Production and factor productivity
    q1      = fill(135.7/12,12),  # "Initial world output 2020 (trill 2019 USD)" 
    al1     = fill(5.84,12),      # "Initial level of total factor productivity"    
    ga1     = fill(0.066,12),     # "Initial growth rate for TFP per 5 years"   
    dela    = fill(0.0015,12),    # "Decline rate of TFP per 5 years"
    gamma   = fill(0.3,12),       # "Capital elasticity in production function"
    # Emission intensity
    gsigma1   = fill(-0.015,12),   # "Initial growth of sigma (per year)"
    delgsig   = fill(0.96,12),     # "Decline rate of gsigma per period"
    asymgsig  = fill(-0.005,12),   # "Asymptotic sigma"
    e1        = fill(37.56/12,12), # "Industrial emissions 2020 (GtCO2 per year)" 
    # Emission from land use
    eland0         = fill(5.9/12,12), # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = fill(0.1,12),    # "Decline rate of land emissions (per period)"  

    # Cost of mitigation
    expcost2  = fill(2.6,12),     # "Exponent of control cost function"

    # Climate change damage
    a2base = fill(0.003467,12), # "Damage quadratic term"

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
    k0       = fill(295/12,12), # "Initial capital stock (10^12 2019 USD)"  
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
    pars = DICEParameters(;regions, weights, gamma, pop1, popadj, popasym, dk, q1, al1, ga1, dela, gsigma1, delgsig, asymgsig, e1, miu1, a1, a2base, a3, expcost2, pback2050, limmiu2070, limmiu2120, limmiu2200, limmiu2300, delmiumax, betaclim, elasmu, prstp, pi_val, k0, siggc1, srf, eland0, deland, f_ghgabate2020, eco2eghgb2020, eco2eghgb2100, emissrat2020, emissrat2100, kwargs...) 
    return pars
end
=#