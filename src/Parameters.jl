"""
In this file we define the struct hosting the parameters of the model and different constructors that define each a "default" of these parameters (e.g. DICE2023 or RICE2023). User can then fine-tune this default with individual parameter choices.
"""




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
    gama    = [0.3]    
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
    optlrsav  = @. (dk + 0.004)/(dk + 0.004*elasmu + rartp)*gama

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
# Constructors

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
  RICE2023()

Parameters constructor with defaults aligned to RICE2023. 

Create a parameters struct with the defaults of the 12-regions RICE2023 model. Different parameters, either the "raw" ones or the "computed ones", can be specified using keywork arguments.

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
    weights = [0.5,0.5],
    gama    = [0.3,0.3], # "Capital elasticity in production function"    
    pop1    = [7752.9/2,7752.9/2 ], # "Initial world population 2020 (millions)" 
    popadj  = [0.145,0.145], # "Growth rate to calibrate to 2050 population projection"
    popasym = [10825/2,10825/2], # "Asymptotic population (millions)"
    dk      = [0.1,0.1], # "Depreciation rate on capital (per year)"     
    q1      = [135.7/2,135.7/2], # "Initial world output 2020 (trill 2019 USD)"   
    al1     = [5.84,5.84], # "Initial level of total factor productivity"    
    ga1     = [0.066,0.066], # "Initial growth rate for TFP per 5 years"   
    dela    = [0.0015,0.0015], # "Decline rate of TFP per 5 years"

    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    gsigma1   = [-0.015,-0.015], # "Initial growth of sigma (per year)"
    delgsig   = [0.96,0.96], # "Decline rate of gsigma per period"
    asymgsig  = [-0.005,-0.005], # "Asymptotic sigma"
    e1        = [37.56/2,37.56/2], # "Industrial emissions 2020 (GtCO2 per year)" 
    miu1      = [0.05,0.05], # "Emissions control rate historical 2020" 
    # --------------------------------------------------------------------
    # Climate damage parameters    
    a1     = [0,0], # "Damage intercept"       
    a2base = [0.003467,0.003467], # "Damage quadratic term"
    a3     = [2,2], # "Damage exponent"
    # --------------------------------------------------------------------
    # Abatement cost
    expcost2  = [2.6,2.6], # "Exponent of control cost function"
    pback2050 = [515,515], # "Cost of backstop in 2019\$ per tCO2 (2050)"
    # --------------------------------------------------------------------
    # Limits on emissions controls
    limmiu2070        = [1.0,1.0], # "Emission control limit from 2070"
    limmiu2120        = [1.1,1.1], # "Emission control limit from 2120"
    limmiu2200        = [1.05,1.05], # "Emission control limit from 2220"
    limmiu2300        = [1.0,1.0], # "Emission control limit from 2300"
    delmiumax         = [0.12,0.12], # "Emission control delta limit per period"
    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing
    betaclim = [0.5,0.5], # "Climate beta"  
    elasmu   = [0.95,0.95], # "Elasticity of marginal utility of consumption" 
    prstp    = [0.001,0.001], # "Pure rate of social time preference"
    pi_val   = [0.05,0.05], # "Capital risk premium (renamed to avoid conflict with Julia's pi)" 
    k0       = [295/2,295/2], # "Initial capital stock (10^12 2019 USD)"  
    siggc1   = [0.01,0.01], # "Annual standard deviation of consumption growth"
    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    srf    = [1000000,1000000], # "Scaling factor for discounting"
    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions  
    eland0         = [5.9/2,5.9/2], # "Carbon emissions from land 2015 (GtCO2 per year)"    
    deland         = [0.1,0.1], # "Decline rate of land emissions (per period)"    
    f_ghgabate2020 = [0.518/2,0.518/2], # "Forcings of abatable non-CO2 GHG in 2020"  
    eco2eghgb2020  = [9.96/2,9.96/2], # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2020"     
    eco2eghgb2100  = [15.5/2,15.5/2], # "Emissions of abatable non-CO2 GHG (GtCO2e) in 2100"   
    emissrat2020   = [1.4,1.4], # "Ratio of CO2e to industrial CO2 in 2020"    
    emissrat2100   = [1.21,1.21], # "Ratio of CO2e to industrial CO2 in 2100"
    kwargs...
    )
    pars = DICEParameters(;regions=regions, weights=weights, gama=gama, pop1=pop1, popadj=popadj, popasym=popasym, dk=dk, q1=q1, al1=al1, ga1=ga1, dela=dela, gsigma1=gsigma1, delgsig=delgsig, asymgsig=asymgsig, e1=e1, miu1=miu1, a1=a1, a2base=a2base, a3=a3, expcost2=expcost2, pback2050=pback2050, limmiu2070=limmiu2070, limmiu2120=limmiu2120, limmiu2200=limmiu2200, limmiu2300=limmiu2300, delmiumax=delmiumax, betaclim=betaclim, elasmu=elasmu, prstp=prstp, pi_val=pi_val, k0=k0, siggc1=siggc1, srf=srf, eland0=eland0, deland=deland, f_ghgabate2020=f_ghgabate2020, eco2eghgb2020=eco2eghgb2020, eco2eghgb2100=eco2eghgb2100, emissrat2020=emissrat2020, emissrat2100=emissrat2100,kwargs...) 
    return pars
end


"""
    @fields_to_vars(t,x)

Utility macro to convert struct fields to local variables (for readibility, so that we can write `parx` instead of using everywhere `p.parx`).
"""
macro fields_to_vars(t::Symbol, x)
    type = Core.eval(__module__, t)
        if !isstructtype(type)
            throw(ArgumentError("@fieldvars only takes struct types, not $type."))
        end 
    esc(:( (; $(fieldnames(type)...)) = $x::$type ))
end