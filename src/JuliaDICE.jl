"""
   JuliaDICE

Implementation of the DICE 2023 model

# Notes:
- Based on DICE2023-b-4-3-10.gms and included files (Nonco2-b-4-3-1.gms and FAIR-beta-4-3-1.gms)
- Variable casing has been harmonized that all parameters and post-optimization computation have lower cases, and all optimization variables have upper case.
"""
module JuliaDICE




export run_dice
using PrecompileTools # just for precompilation
using JuMP, Ipopt


######################################################################
# Raw (exogenous) parameters
######################################################################

# --------------------------------------------------------------------
#  
# Computation of the default miuup (limits on emissions controls)
# Note: `miuup` can be overriden by the user in the `run_dice` function

limmiu2070 = 1.0    # Emission control limit from 2070
limmiu2120 = 1.1  # Emission control limit from 2120
limmiu2200 = 1.05 # Emission control limit from 2220
limmiu2300 = 1.0    # Emission control limit from 2300
delmiumax  = 0.12 # Emission control delta limit per period
tidx  = 1:81
miuupdef       = fill(1.0, 81)
miuupdef[1]       = 0.05; miuupdef[2] = 0.10
@. miuupdef[3:8]  =  delmiumax*(tidx[3:8] - 1)
@. miuupdef[9:11] =  0.85+.05*(tidx[9:11]-8)
miuupdef[12:20]  .= limmiu2070
miuupdef[21:37]  .= limmiu2120
miuupdef[38:57]  .= limmiu2200
miuupdef[58:end] .= limmiu2300


"""
   RawParameters

Raw (exogenous) parameters.
This structure hosts the "default" parameters that can then be modified using keyword arguments in the `run_dice(pars)` function.
"""
@kwdef mutable struct RawParameters
    tstep::Int64   = 5  # Years per period
    ntsteps::Int64 = 81 # Number of time periods

    # --------------------------------------------------------------------
    # Population and technology

    gama::Float64    = 0.3     # Capital elasticity in production function
    pop1::Float64    = 7752.9  # Initial world population 2020 (millions)
    popadj::Float64  = 0.145   # Growth rate to calibrate to 2050 population projection
    popasym::Float64 = 10825   # Asymptotic population (millions)
    dk::Float64      = 0.1     # Depreciation rate on capital (per year)
    q1::Float64      = 135.7   # Initial world output 2020 (trill 2019 USD)
    al1::Float64     = 5.84    # Initial level of total factor productivity
    ga1::Float64     = 0.066   # Initial growth rate for TFP per 5 years
    dela::Float64    = 0.0015  # Decline rate of TFP per 5 years

    # --------------------------------------------------------------------
    # Emissions parameters and Non-CO2 GHG with sigma = emissions/output

    gsigma1::Float64   = -0.015 # Initial growth of sigma (per year)
    delgsig::Float64   = 0.96   # Decline rate of gsigma per period
    asymgsig::Float64  = -0.005 # Asymptotic gsigma
    e1::Float64        = 37.56  # Industrial emissions 2020 (GtCO2 per year)
    miu1::Float64      = 0.05   # Emissions control rate historical 2020
    fosslim::Float64   = 6000   # Maximum cumulative extraction fossil fuels (GtC)
    cumemiss0::Float64 = 633.5  # Cumulative emissions 2020 (GtC)

    # --------------------------------------------------------------------
    # Climate damage parameters

    a1::Float64     = 0        # Damage intercept
    a2base::Float64 = 0.003467 # Damage quadratic term
    a3::Float64     = 2        # Damage exponent

    # --------------------------------------------------------------------
    # Abatement cost

    expcost2::Float64  = 2.6    # Exponent of control cost function
    pback2050::Float64 = 515    # Cost of backstop in 2019$ per tCO2 (2050)
    gback::Float64     = -0.012 # Initial cost decline of backstop cost per year
    cprice1::Float64   = 6      # Carbon price in 2020 (2019$ per tCO2)
    gcprice::Float64   = 0.025  # Growth rate of base carbon price per year

    # --------------------------------------------------------------------
    # Limits on emissions controls
    miuup::Vector{Float64} = miuupdef # Upper bounds on miu

    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing

    betaclim::Float64 = 0.5   # Climate beta
    elasmu::Float64   = 0.95  # Elasticity of marginal utility of consumption
    prstp::Float64    = 0.001 # Pure rate of social time preference
    pi_val::Float64   = 0.05  # Capital risk premium (renamed to avoid conflict with Julia's pi)
    k0::Float64       = 295   # Initial capital stock (10^12 2019 USD)
    siggc1::Float64   = 0.01  # Annual standard deviation of consumption growth

    # --------------------------------------------------------------------
    # Scaling so that MU(C(1)) = 1 and objective function = PV consumption

    srf::Float64    = 1000000    # Scaling factor for discounting
    scale1::Float64 = 0.00891061 # Multiplicative scaling coefficient
    scale2::Float64 = -6275.91   # Additive scaling coefficient

    # --------------------------------------------------------------------
    # Parameters for non-industrial emissions

    eland0::Float64         = 5.9     # Carbon emissions from land 2015 (GtCO2 per year)
    deland::Float64         = 0.1     # Decline rate of land emissions (per period)
    f_misc2020::Float64     = -0.054  # Non-abatable forcings 2020
    f_misc2100::Float64     = 0.265   # Non-abatable forcings 2100
    f_ghgabate2020::Float64 = 0.518   # Forcings of abatable non-CO2 GHG
    f_ghgabate2100::Float64 = 0.957   # Forcings of abatable non-CO2 GHG

    eco2eghgb2020::Float64  = 9.96    # Emissions of abatable non-CO2 GHG (GtCO2e) in 2020
    eco2eghgb2100::Float64  = 15.5    # Emissions of abatable non-CO2 GHG (GtCO2e) in 2100
    emissrat2020::Float64   = 1.4     # Ratio of CO2e to industrial CO2 in 2020
    emissrat2100::Float64   = 1.21    # Ratio of CO2e to industrial CO2 in 2100
    fcoef1::Float64         = 0.00955 # Coefficient of non-CO2 abateable emissions
    fcoef2::Float64         = 0.861   # Coefficient of non-CO2 abateable emissions

    # --------------------------------------------------------------------
    # Parameters for the DFAIR model

    yr0::Int64      = 2020        # Calendar year that corresponds to model year zero

    emshare0::Float64 = 0.2173    # Carbon emissions share into Reservoir 0
    emshare1::Float64 = 0.224     # Carbon emissions share into Reservoir 1
    emshare2::Float64 = 0.2824    # Carbon emissions share into Reservoir 2
    emshare3::Float64 = 0.2763    # Carbon emissions share into Reservoir 3
    tau0::Float64     = 1000000   # Decay time constant for Reservoir 0
    tau1::Float64     = 394.4     # Decay time constant for Reservoir 1
    tau2::Float64     = 36.53     # Decay time constant for Reservoir 2
    tau3::Float64     = 4.304     # Decay time constant for Reservoir 3

    teq1::Float64     = 0.324     # Thermal equilibration parameter for box 1
    teq2::Float64     = 0.44      # Thermal equilibration parameter for box 2
    d1::Float64       = 236       # Thermal response timescale for deep ocean
    d2::Float64       = 4.07      # Thermal response timescale for upper ocean

    irf0::Float64     = 32.4      # Pre-industrial IRF100
    irc::Float64      = 0.019     # Increase in IRF100 with cumulative carbon uptake
    irt::Float64      = 4.165     # Increase in IRF100 with warming
    fco22x::Float64   = 3.93      # Forcings of equilibrium CO2 doubling

    # --------------------------------------------------------------------
    # Initial conditions to be calibrated to history calibration

    mat0::Float64   = 886.5128014 # Initial concentration in atmosphere in 2020 (GtC)

    res00::Float64  = 150.093     # Initial concentration in Reservoir 0 in 2020 (GtC)
    res10::Float64  = 102.698     # Initial concentration in Reservoir 1 in 2020 (GtC)
    res20::Float64  = 39.534      # Initial concentration in Reservoir 2 in 2020 (GtC)
    res30::Float64  = 6.1865      # Initial concentration in Reservoir 3 in 2020 (GtC)

    mateq::Float64  = 588         # Equilibrium concentration in atmosphere (GtC)
    tbox10::Float64 = 0.1477      # Initial temperature box 1 change in 2020 (°C)
    tbox20::Float64 = 1.099454    # Initial temperature box 2 change in 2020 (°C)
    # Changed to prevent numerical instability (GAMS original: 1.24715)...
    tatm0::Float64  = 1.247154    # Initial atmospheric temperature change in 2020 (°C)

end

"""
    @fields_to_vars(t,x)

Utility macro to convert a struct fields to local variables (for readibility, so that we can write `parx` instead of using everywhere `p.parx`).
"""
macro fields_to_vars(t::Symbol, x)
    type = Core.eval(__module__, t)
        if !isstructtype(type)
            throw(ArgumentError("@fieldvars only takes struct types, not $type."))
        end 
    esc(:( (; $(fieldnames(type)...)) = $x::$type ))
end

"""
    run_dice(;optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5), bounds=Dict{String,Tuple{String,String}}(),kwargs...)

Run the DICE model (currently v 2023).

This function runs the DICE model and returns the results as a named tuple.

# Parameters

- `optimizer`: The optimizer to use and eventually its options. Defaults: [`optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5)`].
- `bounds`: A dictionary of equality or inequality constraints. Each constraint should be specified with the variable name as the key, and a two-elements tuple as the value. The first element is either "<=", ">=", or "==", and the second element is the right-hand side of the constraint (single value or a vector of ntimesteps length). Default: (empty dictionary). See the source code for the name of the variables of the model.
- `kwargs`: Keyword arguments to override the default parameter values. Again, see the model source code for the name of the parameters.

# Outputs

- A named tuple containing the following fields: `solved`, `status`, `times`, `tidx`, the post_process computed values and the optimization variables.

# Examples:

```Julia
res = run_dice()
ECO2_opt = res.ECO2
plot(res.times[1:11] .+ 2020,ECO2_opt[1:11],ylim=(0,80), title="CO₂ emissions",ylabel="GtCO₂/yr",label="C/B optimal", markershape=:circle, markercolor=:white)
```

```Julia
res_crazy = run_dice(optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0), bounds = Dict("MIU"=>("==",1.0), "TATM"=>("<=",15), "Y" =>(">=",[fill(floatmin(Float64),10);fill(0.1,71)]), "ECO2" =>("<=",10000)), a2base = 0.01)
```


# Notes
- The `bounds` adds constraint to the problem, but do not substitute to hard-written bounds in the model. In particular, for the upper limit of the emissions controls, use the parameter `miuup` instead 



"""
function run_dice(;optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0), bounds=Dict{String,Tuple{String,String}}(),kwargs...) 

    issubset(keys(kwargs), fieldnames(RawParameters)) || error("Not all keywords are valid parameters.")
    p = RawParameters(;kwargs...)   # Override the default parameter values with the keyword arguments
    @fields_to_vars RawParameters p # Copy of the RowParameters fields to local variables (for readibility)

    ######################################################################
    # Computed parameters
    ######################################################################

    # --------------------------------------------------------------------
    # Preferences, growth uncertainty, and timing

    # Time periods sequence
    times = 0:tstep:(ntsteps*tstep)-1  # 0,5,10,...,400
    # Time periods index sequences
    tidx  = 1:ntsteps                 # 1,2,3,...,81
    t0idx = 0:ntsteps-1               # 0,1,2,...,80

    # Risk-adjusted rate of time preference
    rartp = exp(prstp + betaclim*pi_val)-1
    
    # --------------------------------------------------------------------
    # Precautionary dynamic parameters

    # Variance of per capita consumption 
    varpcc           = [min((siggc1^2)*t,(siggc1^2)*tstep*47) for t in times]
    # Precautionary rate of return   
    rprecaut         = @. -0.5 * varpcc * elasmu^2
    # STP factor without precautionary factor
    rr1              = @. 1 / ((1+rartp)^times)
    # STP with precautionary factor
    rr               = @. rr1 * (1+rprecaut) ^ -times
    # Level of population and labor
    l                = fill(pop1, ntsteps); [l[ti] = l[ti-1] * (popasym / l[ti-1])^popadj for ti in tidx[2:end]]
    # Growth rate of Total Factor Productivity
    ga               = @. ga1*exp(-dela*times)
    # Level of total factor productivity
    al               = fill(al1, ntsteps); [al[ti] = al[ti-1] / (1 - ga[ti-1]) for ti in tidx[2:end]]
    # Optimal long-run savings rate used for transversality
    optlrsav         = (dk + 0.004)/(dk + 0.004*elasmu + rartp)*gama
    # Carbon price in base case
    cpricebase       = @. cprice1*(1+gcprice)^times
    # Backstop price 2019$ per ton CO2
    pbacktime        = @. pback2050*exp(-tstep*0.01*(tidx-7))
    pbacktime[8:end] = @. pback2050*exp(-tstep*0.001*(tidx[8:end]-7))
    # Carbon intensity 2020 kgCO2-output 2020
    sig1             = e1/(q1*(1-miu1))

    # Change in sigma (rate of decarbonization)
    gsig             = @. min(gsigma1*delgsig ^(tidx-1),asymgsig)

    # CO2-emissions output ratio
    sigma            = fill(sig1,ntsteps)
    [sigma[ti] = sigma[ti-1] * exp(tstep*gsig[ti-1]) for ti in tidx[2:end] ]

     # ------------------------------------------------------------------------------
    # Parameters emissions and non-CO2 
    
    eland          = @. eland0*(1-deland)^t0idx  
    co2e_ghgabateb = eco2eghgb2020 .+ [(eco2eghgb2100-eco2eghgb2020) * min(1,ti/16) for ti in t0idx]
    f_misc         = f_misc2020    .+ [(f_misc2100-f_misc2020) * min(1,ti/16) for ti in t0idx]
    emissrat       = emissrat2020  .+ [(emissrat2100-emissrat2020) * min(1,ti/16) for ti in t0idx]
    sigmatot       = @. sigma * emissrat
    cost1tot       = @. pbacktime * sigmatot / expcost2/1000

    ######################################################################
    # Optimization model & computation options
    ######################################################################


    m = Model(optimizer)


    ######################################################################
    # Variables declaration
    ######################################################################

    @variables m begin

        ECO2[tidx]             # Total CO2 emissions (GtCO2 per year)
        ECO2E[tidx]            # Total CO2e emissions including abateable nonCO2 GHG (GtCO2 per year)
        EIND[tidx]             # Industrial CO2 emissions (GtCO2 per yr)
        # TODO: check this, as in GAMS it is declared without time index, but then it is used with time index!
        F_GHGABATE[tidx]       # Forcings abateable nonCO2 GHG    
        
        MIU[tidx] >= 0.0       # Emission control rate GHGs (**control**)
        C[tidx] >= 0.0         # Consumption (trillions 2019 US dollars per year)
        K[tidx] >= 0.0         # Capital stock (trillions 2019 US dollars)
        CPC[tidx]              # Per capita consumption (thousands 2019 USD per year)
        I[tidx] >= 0.0         # Investment (trillions 2019 USD per year)
        S[tidx]                # Gross savings rate as fraction of gross world product (**control**)
        Y[tidx] >= 0.0         # Gross world product net of abatement and damages (trillions 2019 USD per year)
        YGROSS[tidx] >= 0.0    # Gross world product GROSS of abatement and damages (trillions 2019 USD per year)
        YNET[tidx] >= 0.0      # Output net of damages equation (trillions 2019 USD per year)
        DAMAGES[tidx]          # Damages (trillions 2019 USD per year)
        DAMFRAC[tidx]          # Damages as fraction of gross output
        ABATECOST[tidx]        # Cost of emissions reductions  (trillions 2019 USD per year)
        MCABATE[tidx]          # Marginal cost of abatement (2019$ per ton CO2)
        CCATOT[tidx]           # Total carbon emissions (GtC)
        PERIODU[tidx]          # One period utility function
        CPRICE[tidx]           # Carbon price (2019$ per ton of CO2)
        TOTPERIODU[tidx]       # Period utility
        UTILITY                # Welfare function

        # --------------------------------------------------------------------
        # Climate variables
        # Note: Stock variables correspond to levels at the END of the period

        FORC[tidx]             # Increase in radiative forcing (watts per m2 from 1765)
        TATM[tidx] >= 0.5      # Increase temperature of atmosphere (degrees C from 1765)     
        TBOX1[tidx]            # Increase temperature of box 1 (degrees C from 1765)
        TBOX2[tidx]            # Increase temperature of box 2 (degrees C from 1765)
        RES0[tidx]             # Carbon concentration in Reservoir 0 (GtC from 1765)
        RES1[tidx]             # Carbon concentration in Reservoir 1 (GtC from 1765)
        RES2[tidx]             # Carbon concentration in Reservoir 2 (GtC from 1765)
        RES3[tidx]             # Carbon concentration in Reservoir 3 (GtC from 1765)
        MAT[tidx] >= 0.0       # Carbon concentration increase in atmosphere (GtC from 1765)
        CACC[tidx]             # Accumulated carbon in ocean and other sinks (GtC)
        IRFT[tidx] >= 0.0      # IRF100 at time t
        ALPHA[tidx] >= 0.0     # Carbon decay time scaling factor (**control**)
    end

    ######################################################################
    # Constraints
    ######################################################################

    # --------------------------------------------------------------------
    # Climate variables
    # Note: initial conditionsand stability bounds are set here as contraints, we can also try as fixed value

    # Reservoir 0 law of motion
    @constraint(m, res0lom[ti in tidx], RES0[ti] ==  ((ti == 1) ? res00 : (emshare0*tau0*ALPHA[ti]*(ECO2[ti]/3.667))*(1-exp(-tstep/(tau0*ALPHA[ti])))+RES0[ti-1]*exp(-tstep/(tau0*ALPHA[ti]))))

    # Reservoir 1 law of motion
    @constraint(m, res1lom[ti in tidx], RES1[ti] == ((ti == 1) ? res10 :  (emshare1*tau1*ALPHA[ti]*(ECO2[ti]/3.667))*(1-exp(-tstep/(tau1*ALPHA[ti])))+RES1[ti-1]*exp(-tstep/(tau1*ALPHA[ti]))))

    # Reservoir 2 law of motion
    @constraint(m, res2lom[ti in tidx], RES2[ti] ==  ((ti == 1) ? res20 : (emshare2*tau2*ALPHA[ti]*(ECO2[ti]/3.667))*(1-exp(-tstep/(tau2*ALPHA[ti])))+RES2[ti-1]*exp(-tstep/(tau2*ALPHA[ti]))))

    # Reservoir 3 law of motion 
    @constraint(m, res3lom[ti in tidx], RES3[ti] ==  ((ti == 1) ? res30 : (emshare3*tau3*ALPHA[ti]*(ECO2[ti]/3.667))*(1-exp(-tstep/(tau3*ALPHA[ti])))+RES3[ti-1]*exp(-tstep/(tau3*ALPHA[ti]))))

    # Atmospheric concentration equation
    @constraint(m, matlb[ti in tidx], MAT[ti] >= 10.0)
    @constraint(m, mmat[ti in tidx], MAT[ti] == ((ti == 1) ? mat0 : mateq + RES0[ti] + RES1[ti] + RES2[ti] + RES3[ti]))

    # Accumulated carbon in sinks equation
    @constraint(m, cacceq[ti in tidx], CACC[ti] == CCATOT[ti]-(MAT[ti]-mateq))

    # Radiative forcing equation
    @constraint(m, force[ti in tidx], FORC[ti] == fco22x*((log((MAT[ti]/mateq))/log(2))) + f_misc[ti]+F_GHGABATE[ti] )

    # Temperature box 1 law of motion
    @constraint(m, tbox1eq[ti in tidx], TBOX1[ti] ==  ((ti == 1) ? tbox10 : TBOX1[ti-1]*exp(-tstep/d1)+teq1*FORC[ti]*(1-exp(-tstep/d1))))

    # Temperature box 2 law of motion
    @constraint(m, tbox2eq[ti in tidx], TBOX2[ti] ==  ((ti == 1) ? tbox20 : TBOX2[ti-1]*exp(-tstep/d2)+teq2*FORC[ti]*(1-exp(-tstep/d2))))

    # Temperature-climate equation for atmosphere
    @constraint(m, tatmlb[ti in tidx], TATM[ti] >= 0.5)
    @constraint(m, tatmub[ti in tidx], TATM[ti] <= 20.0)
    @constraint(m, tatmeq[ti in tidx], TATM[ti] ==  ((ti == 1) ? tatm0 : TBOX1[ti]+TBOX2[ti]))

    # Left-hand side of IRF100 equation
    @constraint(m, irfeqlhs[ti in tidx],  IRFT[ti] == (ALPHA[ti]*emshare0*tau0*(1-exp(-100/(ALPHA[ti]*tau0))))+(ALPHA[ti]*emshare1*tau1*(1-exp(-100/(ALPHA[ti]*tau1))))+(ALPHA[ti]*emshare2*tau2*(1-exp(-100/(ALPHA[ti]*tau2))))+(ALPHA[ti]*emshare3*tau3*(1-exp(-100/(ALPHA[ti]*tau3)))))

    # Right-hand side of IRF100 equation
    @constraint(m, irfeqrhs[ti in tidx],  IRFT[ti] == irf0+irc*CACC[ti]+irt*TATM[ti])

    # --------------------------------------------------------------------
    # Emissions and Damages

    # CO2 Emissions equation
    @constraint(m, eco2eq[ti in tidx], ECO2[ti] == (sigma[ti]*YGROSS[ti] + eland[ti])*(1-MIU[ti]))

    # Industrial CO2 equation
    @constraint(m, eindeq[ti in tidx], EIND[ti] == (sigma[ti]*YGROSS[ti])*(1-MIU[ti]))

    # CO2E Emissions equation
    @constraint(m, eco2eeq[ti in tidx], ECO2E[ti] == (sigma[ti]*YGROSS[ti] + eland[ti] + co2e_ghgabateb[ti])*(1-MIU[ti]))

    # Forcings abateable nonCO2 GHG equation 
    @constraint(m, f_ghgabateeq[ti in tidx], F_GHGABATE[ti] == ((ti == 1) ? f_ghgabate2020 : fcoef2*F_GHGABATE[ti-1]+ fcoef1*co2e_ghgabateb[ti-1]*(1-MIU[ti-1])))

    # --------------------------------------------------------------------
    # Emissions and damage

    # Cumulative total carbon emissions
    @constraint(m, ccatoteq[ti in tidx], CCATOT[ti] == ((ti==1) ? cumemiss0 : CCATOT[ti-1] +  ECO2[ti-1]*(tstep/3.666)))

    # Equation for damage fraction
    @constraint(m, damfraceq[ti in tidx], DAMFRAC[ti] == (a1*TATM[ti])+(a2base*TATM[ti]^a3))

    # Damage equation
    @constraint(m, dameq[ti in tidx], DAMAGES[ti] == YGROSS[ti] * DAMFRAC[ti])

    # Cost of emissions reductions equation
    @constraint(m, abateeq[ti in tidx], ABATECOST[ti] == YGROSS[ti] * cost1tot[ti] * (MIU[ti]^expcost2))

    # Equation for MC abatement
    @constraint(m, mcabateeq[ti in tidx], MCABATE[ti] == pbacktime[ti] * MIU[ti]^(expcost2-1))

    # Carbon price equation from abatement
    @constraint(m, carbpriceeq[ti in tidx], CPRICE[ti] == pbacktime[ti] * MIU[ti]^(expcost2-1))

    # --------------------------------------------------------------------
    # Economic variables

    # Output gross equation
    @constraint(m, ygrosseq[ti in tidx], YGROSS[ti] == (al[ti]*(l[ti]/1000)^(1-gama))*(K[ti]^gama))

    # Output net of damage equation
    @constraint(m, yneteq[ti in tidx], YNET[ti] == YGROSS[ti]*(1-DAMFRAC[ti]))

    # Output net equation
    @constraint(m, yy[ti in tidx], Y[ti] == YNET[ti] - ABATECOST[ti])

    # Consumption equation
    @constraint(m, cc[ti in tidx], C[ti] == Y[ti] - I[ti])
    @constraint(m, clb[ti in tidx],  C[ti] >= 2)

    # Per capita consumption definition
    @constraint(m, cpce[ti in tidx], CPC[ti] == 1000 * C[ti] / l[ti])
    @constraint(m, cpclb[ti in tidx],  CPC[ti] >= 0.01)

    # Savings rate equation
    @constraint(m, seq[ti in tidx], I[ti] == S[ti] * Y[ti])

    # Capital balance equation
    @constraint(m, kk0, K[1] == k0)
    @constraint(m, kk[ti in tidx[2:end]],  K[ti] <= (1-dk)^tstep * K[ti-1] + tstep * I[ti-1])
    @constraint(m, klb[ti in tidx],  K[ti] >= 1)

    # -------------------------------------------------------------------- 
    # Utility

    # Instantaneous utility function equation
    @constraint(m, periodueq[ti in tidx], PERIODU[ti] == ((C[ti]*1000/l[ti])^(1-elasmu)-1)/(1-elasmu)-1)

    # Period utility
    @constraint(m, totperiodueq[ti in tidx], TOTPERIODU[ti] == PERIODU[ti] * l[ti] * rr[ti])

    # Total utility
    @constraint(m, utileq, UTILITY == tstep * scale1 * sum(TOTPERIODU[ti] for ti in tidx) + scale2)


    # -------------------------------------------------------------------- 
    # Other upper and lower bounds for stability

    @constraint(m, alphalb[ti in tidx] , ALPHA[ti] >= 0.1)
    @constraint(m, alphaub[ti in tidx] , ALPHA[ti] <= 100.0)
    @constraint(m, miuub[ti in tidx] , MIU[ti] <= miuup[ti])
    @constraint(m, sfix[ti in tidx[38:end]] , S[ti] == 0.28)

    # -------------------------------------------------------------------- 
    # Scenario-dependant constraints

    mvars  = object_dictionary(m)
    for (k,v) in bounds
        v_vector = (ndims(v[2]) == 0) ? fill(v[2],ntsteps) : v[2]
        if (v[1] == "<=")
            scenboundseq = @constraint(m, [ti in tidx], mvars[Symbol(k)][ti] <= v_vector[ti])
        elseif (v[1] == ">=")
            scenboundseq = @constraint(m, [ti in tidx], mvars[Symbol(k)][ti] >= v_vector[ti])
        elseif (v[1] == "==")
            scenboundseq = @constraint(m, [ti in tidx], mvars[Symbol(k)][ti] == v_vector[ti])
        else
            error("Unknown constraint type $(v[1])")
        end
    end

    #####################################################################
    # OBJECTIVE FUNCTION
    #####################################################################
    @objective(m, Max, UTILITY)

    #####################################################################
    # OPTIMIZATION
    #####################################################################
    optimize!(m)

    status = termination_status(m)

    if(!is_solved_and_feasible(m))
        @warn "Model not solved"
        return((solved=false,status=status))
    end

    #####################################################################
    # POST-OPTIMIZATION COMPUTATIONS
    #####################################################################

    rfactlong  =  collect(@. srf *(value(CPC)/value(CPC[1]))^(-elasmu)*rr)
    rlong      =  [-log(rfactlong[ti]/srf)/(tstep*ti) for ti in tidx]
    rshort     =  vcat(missing,[-log(rfactlong[ti]/rfactlong[ti-1])/tstep for ti in tidx[2:end]])

    scc        = collect(@. -1000*dual(eco2eq) /(0.00001+dual(cc)))
    ppm        = collect(value.(MAT) ./ 2.13)
    abaterat   = collect(@. value(ABATECOST)/value(Y))
    atfrac2020 = collect(@. (value(MAT)-mat0)/(value(CCATOT)+0.00001-cumemiss0))
    atfrac1765 = collect(@. (value(MAT)-mateq)/(value(CCATOT)+0.00001))
    forc_co2   = collect(@. fco22x*log(value(MAT)/mateq)/log(2))

    # Return results as named tuple
    return (solved=true, status=status, times=times, tidx=tidx, rlong=rlong, rshort=rshort, scc=scc, ppm=ppm, abaterat=abaterat, atfrac2020=atfrac2020, atfrac1765=atfrac1765, forc_co2=forc_co2, ECO2=collect(value.(ECO2)), ECO2E=collect(value.(ECO2E)), EIND=collect(value.(EIND)), F_GHGABATE=collect(value.(F_GHGABATE)),MIU=collect(value.(MIU)), C=collect(value.(C)), K=collect(value.(K)), CPC=collect(value.(CPC)), I=collect(value.(I)), S=collect(value.(S)), Y=collect(value.(Y)), YGROSS=collect(value.(YGROSS)), YNET=collect(value.(YNET)), DAMAGES=collect(value.(DAMAGES)), DAMFRAC=collect(value.(DAMFRAC)), ABATECOST=collect(value.(ABATECOST)), MCABATE=collect(value.(MCABATE)), CCATOT=collect(value.(CCATOT)), PERIODU=collect(value.(PERIODU)), CPRICE=collect(value.(CPRICE)), TOTPERIODU=collect(value.(TOTPERIODU)), UTILITY=value(UTILITY), FORC=collect(value.(FORC)), TATM=collect(value.(TATM)), TBOX1=collect(value.(TBOX1)), TBOX2=collect(value.(TBOX2)), RES0=collect(value.(RES0)), RES1=collect(value.(RES1)), RES3=collect(value.(RES3)), MAT=collect(value.(MAT)), CACC=collect(value.(CACC)), IRFT=collect(value.(IRFT)), ALPHA=collect(value.(ALPHA)))
end


include("Precompilation.jl") 

end # module JuliaDICE
