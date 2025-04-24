"""
Part of [DICEModel](https://github.com/sylvaticus/DICEModel.jl). Licence is MIT.

This file contains `run_dice`, the low-level function that runs the optimization
"""

"""
    run_dice(pars;optimizer,bounds)
    run_dice(;optimizer,bounds,kwargs...)

Run the D(R)ICE models (currently with the structure and, by default, the data of DICE2023), possibly with custom optimiser, bounds or parameters.

This function runs the DICE model and returns the results as a named tuple. Note that starting from DICEModel v0.2, a regional dimension is always present, and DICE is simply treated as RICE with a single region.  

# Function arguments

## Positional:
- `pars`: An istance of the [`DICEModel`](@ref) struct containing the needed parameters

## Keyword arguments:
- `optimizer': The optimiser to use and possibly its options. Defaults to: [`optimizer = optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0, "max_wall_time"=>10.0^20, "max_cpu_time" => 10.0^20, "max_iter" => 3000, "acceptable_tol" =>10^-6, "acceptable_iter" => 15, "acceptable_dual_inf_tol" =>10.0^10, "acceptable_constr_viol_tol" => 0.01, "acceptable_compl_inf_tol" =>0.01, "acceptable_obj_change_tol" =>10.0^20)`]. All, except the print levels, are the Ipopt defaults.
- `bounds``: A dictionary of equality or inequality constraints. Each constraint should be specified with the variable name as key and a two-element tuple as value. The first element is either "<=", ">=" or "==", and the second element is the right-hand side of the constraint (a single value, a vector of ntimesteps length or a matrix of ntimesteps x nregions). Default: (empty dictionary). See the [source code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/DICEModel.jl) for the names of the model variables.
- `kwargs``: Keyword arguments to override the default parameter values of DICE2023. See the documentation for the [`DICEParameters`](@ref) structure for the available model parameters. 

**WARNING**: Sometimes changing a parameter doesn't lead to the expected behavior. This is because the model (in its original GAMS form that has been re-implemented in this package) performs some calibrations with the parameters, so several parameters have to be changed together. For example all the scenarios that test different discount rates don't change only the `prstp` parameter, but compute several other parameters, sometimes in a different matter than the default model, and have different calibration for initial conditions. Always check the [source code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/CoreModel.jl) to make sure that the parameter you want to change doesn't have other side effects in the model.

# Outputs

- A named tuple containing the following fields: `solved`, `status`, `times`, `tidx`, the post_process computed values, the optimisation variables, the parameters structure (`pars`).

# Examples:

```Julia
res = run_dice()
ECO2_opt = res.ECO2
plot(res.times[1:11] .+ 2020,ECO2_opt[1:11],ylim=(0,80), title="CO₂ emissions",ylabel="GtCO₂/yr",label="C/B optimal", markershape=:circle, markercolor=:white)
```

```Julia
res_crazy = run_dice(optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0), bounds = Dict("MIU"=>("==",1.0), "TATM"=>("<=",15), "Y" =>(">=",[fill(floatmin(Float64),10);fill(0.1,71)]), "ECO2" =>("<=",10000)), a2base = 0.01)
```

```julia
w_dev_country_priority = [1,1,1,1.5,2,2,3,3,2,5,5,5] # Utility weights by region
res_12regions_devprior = run_dice(RICE2023(;weights=w_dev_country_priority))
```

# Notes
- The `bounds` add constraints to the problem, but do not replace hard written bounds in the model. In particular, the `miuup` parameter should be used instead for the upper limit of emission controls.
- Bounds are always intended for the full time steps. If you need a bound for a subset of time steps (e.g. the first time step), you still need to assemble your full time array of the bound using `floatmin(Float64)` or `floatmax(Float64)` as appropriate.
- The version with keywords arguments is a tiny wrapper (calls) the version with the parameters struct, itself built using the `DICE2023` function with the provided keyword arguments.
- The weights used in the [`DICEParameters`](@ref) struct are exogenous, are NOT the Negishi weights. The `run_dice` function can eventually be used iteractively to look for these weights.
"""
function run_dice(
      pars::DICEParameters;
      optimizer = optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0, "max_wall_time"=>10.0^20, "max_cpu_time" => 10.0^20, "max_iter" => 3000, "acceptable_tol" =>10^-6, "acceptable_iter" => 15, "acceptable_dual_inf_tol" =>10.0^10, "acceptable_constr_viol_tol" => 0.01, "acceptable_compl_inf_tol" =>0.01, "acceptable_obj_change_tol" =>10.0^20),
      bounds = Dict{String,Tuple{String,String}}()
    ) 

    Random.seed!(123)
    
    @fields_to_vars DICEParameters pars # Copy of the RowParameters fields to local variables (for readibility)

    weights = scaleweights(weights)

    ######################################################################
    # Optimization model & computation options
    ######################################################################

    m = Model(optimizer)

    ######################################################################
    # Variables declaration
    ######################################################################

    @variables m begin
        EIND_R[tidx,ridx]       # Industrial CO2 emissions (GtCO2 per yr), regional
        EIND[tidx]              # Industrial CO2 emissions (GtCO2 per yr)
        ECO2_R[tidx,ridx]       # Total CO2 emissions (GtCO2 per year), regional
        ECO2[tidx]              # Total CO2 emissions (GtCO2 per year)
        ECO2E_R[tidx,ridx]      # Total CO2e emissions including abateable nonCO2 GHG (GtCO2 per year), regional
        ECO2E[tidx]             # Total CO2e emissions including abateable nonCO2 GHG (GtCO2 per year)
        F_GHGABATE_R[tidx,ridx] # Forcings abateable nonCO2 GHG, regional   
    
        
        MIU_R[tidx,ridx] >= 0.0      # Emission control rate GHGs (**control**), regional
        C_R[tidx,ridx]  >= 0.0       # Consumption (trillions 2019 US dollars per year), regional
        C[tidx] >= 0.0         # Consumption (trillions 2019 US dollars per year)
        K_R[tidx,ridx] >= 0.0         # Capital stock (trillions 2019 US dollars), regional
        K[tidx] >= 0.0         # Capital stock (trillions 2019 US dollars)
        CPC_R[tidx,ridx]              # Per capita consumption (thousands 2019 USD per year), regional
        CPC[tidx]              # Per capita consumption (thousands 2019 USD per year)
        I_R[tidx,ridx] >= 0.0         # Investment (trillions 2019 USD per year), regional
        I[tidx] >= 0.0         # Investment (trillions 2019 USD per year)
        S_R[tidx,ridx]               # Gross savings rate as fraction of gross world product (**control**), regional
        S[tidx]                # Gross savings rate as fraction of gross world product (**control**)
        Y_R[tidx,ridx] >= 0.0         # Gross world product net of abatement and damages (trillions 2019 USD per year), regional
        Y[tidx] >= 0.0         # Gross world product net of abatement and damages (trillions 2019 USD per year)
        YGROSS_R[tidx,ridx] >= 0.0    # Gross world product GROSS of abatement and damages (trillions 2019 USD per year), regional
        YGROSS[tidx] >= 0.0    # Gross world product GROSS of abatement and damages (trillions 2019 USD per year)
        YNET_R[tidx,ridx] >= 0.0      # Output net of damages equation (trillions 2019 USD per year)
        YNET[tidx] >= 0.0      # Output net of damages equation (trillions 2019 USD per year)
        0.0 <= DAMFRAC_R[tidx,ridx] <= 1.0    # Damages as fraction of gross output, regional
        DAMAGES[tidx]          # Damages (trillions 2019 USD per year)
        DAMAGES_R[tidx,ridx] >= 0.0         # Damages (trillions 2019 USD per year)
        ABATECOST_R[tidx, ridx]        # Cost of emissions reductions  (trillions 2019 USD per year), regional
        ABATECOST[tidx]        # Cost of emissions reductions  (trillions 2019 USD per year)
        CCATOT[tidx]           # Total carbon emissions (GtC)
        PERIODU_R[tidx,ridx]   # One period utility function, regional
        PERIODU[tidx]          # One period utility function
        CPRICE_R[tidx,ridx]           # Carbon price (2019$ per ton of CO2) from marginal abatement cost
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

    #Reservoir 1 law of motion
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
    @constraint(m, force[ti in tidx], FORC[ti] == fco22x*(log(MAT[ti]/mateq)/log(2)) + f_misc[ti]+sum(F_GHGABATE_R[ti,ri] for ri in ridx))

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

    # Industrial CO2 equation, regional
    @constraint(m, eindeq_r[ti in tidx, ri in ridx], EIND_R[ti,ri] == (sigma[ti,ri]*YGROSS_R[ti,ri])*(1-MIU_R[ti,ri]))
    # Industrial CO2 equation, total
    @constraint(m, eindeq[ti in tidx], EIND[ti] == sum(EIND_R[ti,ri] for ri in ridx))


    # CO2 Emissions equation, regional
    @constraint(m, eco2eq_r[ti in tidx, ri in ridx], ECO2_R[ti,ri] == ((sigma[ti,ri]*YGROSS_R[ti,ri]) + eland[ti,ri]) *(1-MIU_R[ti,ri]))
    # CO2 Emissions equation
    @constraint(m, eco2eq[ti in tidx], ECO2[ti] == sum(ECO2_R[ti,ri] for ri in ridx))

    # CO2E Emissions equation, regional
    @constraint(m, eco2eeq_r[ti in tidx, ri in ridx], ECO2E_R[ti,ri] == ((sigma[ti,ri]*YGROSS_R[ti,ri]) + eland[ti,ri] + co2e_ghgabateb[ti,ri])*(1-MIU_R[ti,ri]))
    # CO2E Emissions equation
    @constraint(m, eco2eeq[ti in tidx], ECO2E[ti] == sum(ECO2E_R[ti,ri] for ri in ridx))

    # Forcings abateable nonCO2 GHG equation, regional 
    @constraint(m, f_ghgabateeq_r[ti in tidx, ri in ridx], F_GHGABATE_R[ti,ri] == ((ti == 1) ? f_ghgabate2020[ri] : fcoef2*F_GHGABATE_R[ti-1,ri]+ fcoef1*co2e_ghgabateb[ti-1,ri]*(1-MIU_R[ti-1,ri])))


    # --------------------------------------------------------------------
    # Emissions and damage

    # Cumulative total carbon emissions
    @constraint(m, ccatoteq[ti in tidx], CCATOT[ti] == ((ti==1) ? cumemiss0 : CCATOT[ti-1] +  ECO2[ti-1]*(tstep/3.666)))

    # Equation for damage fraction
    @constraint(m, damfraceq_r[ti in tidx, ri in ridx], DAMFRAC_R[ti,ri] == (a1[ri]*TATM[ti])+(a2base[ri]*TATM[ti]^a3[ri]))

    # Damage equation
    @constraint(m, dameq_r[ti in tidx, ri in ridx], DAMAGES_R[ti,ri] == YGROSS_R[ti,ri] * DAMFRAC_R[ti,ri])
    @constraint(m, dameq[ti in tidx], DAMAGES[ti] == sum(DAMAGES_R[ti,ri] for ri in ridx))

    # Cost of emissions reductions equation
    @constraint(m, abateeq_r[ti in tidx,ri in ridx], ABATECOST_R[ti,ri] == YGROSS_R[ti,ri] * cost1tot[ti,ri] * (MIU_R[ti,ri]^expcost2[ri]))
    @constraint(m, abateeq[ti in tidx], ABATECOST[ti] == sum(ABATECOST_R[ti,ri] for ri in ridx))

    # Carbon price equation from abatement
    @constraint(m, carbpriceeq_r[ti in tidx, ri in ridx], CPRICE_R[ti,ri] == pbacktime[ti,ri] * MIU_R[ti,ri]^(expcost2[ri]-1))

    # --------------------------------------------------------------------
    # Economic variables

    # Output gross equation
    @constraint(m, ygrosseq_r[ti in tidx, ri in ridx], YGROSS_R[ti,ri] == (al[ti,ri]*(l[ti,ri]/1000)^(1-gamma[ri]))*(K[ti]^gamma[ri]))
    @constraint(m, ygrosseq[ti in tidx], YGROSS[ti] ==  sum(YGROSS_R[ti,ri] for ri in ridx))

    # Output net of damage equation
    @constraint(m, yneteq_r[ti in tidx, ri in ridx], YNET_R[ti,ri] == YGROSS_R[ti,ri]*(1-DAMFRAC_R[ti,ri]))
    @constraint(m, yneteq[ti in tidx], YNET[ti] == sum(YNET_R[ti,ri] for ri in ridx))

    # Output net equation
    @constraint(m, yy_r[ti in tidx, ri in ridx], Y_R[ti,ri] == YNET_R[ti,ri] - ABATECOST_R[ti,ri])
    @constraint(m, yy[ti in tidx], Y[ti] == sum(Y_R[ti,ri] for ri in ridx))

    # Consumption equation
    @constraint(m, cc_r[ti in tidx,ri in ridx], C_R[ti,ri] == Y_R[ti,ri] - I_R[ti,ri])
    @constraint(m, cc[ti in tidx], C[ti] == sum(C_R[ti,ri] for ri in ridx))
    

    # Per capita consumption definition
    @constraint(m, cpce_r[ti in tidx, ri in ridx], CPC_R[ti,ri] == 1000 * C_R[ti,ri] / l[ti,ri])
    @constraint(m, cpce[ti in tidx], CPC[ti] == 1000 * C[ti]  / sum(l[ti,ri] for ri in ridx))
    

    # Savings rate equation
    @constraint(m, ieq_r[ti in tidx, ri in ridx], I_R[ti,ri] == S_R[ti,ri] * Y_R[ti,ri])
    @constraint(m, ieq[ti in tidx], I[ti] == sum(I_R[ti,ri] for ri in ridx))
    @constraint(m, seq[ti in tidx], S[ti] == I[ti] / Y[ti]) 

    # Capital balance equation
    @constraint(m, kk0_r[ri in ridx], K_R[1,ri] == k0[ri])
    @constraint(m, kk_r[ti in tidx[2:end],ri in ridx],  K_R[ti,ri] <= (1-dk[ri])^tstep * K_R[ti-1,ri] + tstep * I_R[ti-1,ri])
    @constraint(m, kk[ti in tidx],  K[ti] == sum(K_R[ti,ri] for ri in ridx))
    @constraint(m, klb_r[ti in tidx, ri in ridx],  K_R[ti,ri] >= 1)

    # -------------------------------------------------------------------- 
    # Utility

    # Instantaneous utility function equation
    @constraint(m, periodueq_r[ti in tidx, ri in ridx], PERIODU_R[ti,ri] == ((C_R[ti,ri]*1000/l[ti,ri])^(1-elasmu[ri])-1)/(1-elasmu[ri])-1)
    @constraint(m, periodueq[ti in tidx], PERIODU[ti] == sum(PERIODU_R[ti,ri] for ri in ridx))

    # Period utility
    @constraint(m, totperiodueq[ti in tidx], TOTPERIODU[ti] == sum(PERIODU_R[ti,ri] * l[ti,ri] * rr[ti,ri] * weights[ri] / sum(weights) for ri in ridx)) # note: here rr will be rr[ti,ri] to consider the regional weigths

    # Total utility
    @constraint(m, utileq, UTILITY == tstep * scale1 * sum(TOTPERIODU[ti] for ti in tidx) + scale2)


    # -------------------------------------------------------------------- 
    # Other upper and lower bounds for stability

    @constraint(m, alphalb[ti in tidx] , ALPHA[ti] >= 0.1)
    @constraint(m, alphaub[ti in tidx] , ALPHA[ti] <= 100.0)
    @constraint(m, miuub[ti in tidx, ri in ridx] , MIU_R[ti,ri] <= miuup[ti,ri])
    @constraint(m, sfix[ti in tidx[38:end]] , S[ti] == 0.28)
    @constraint(m, clb[ti in tidx,ri in ridx],  C_R[ti,ri] >= 2.0)
    @constraint(m, cpclb[ti in tidx, ri in ridx],  CPC_R[ti,ri] >= 0.01)

    # -------------------------------------------------------------------- 
    # Scenario-dependant constraints

    mvars  = object_dictionary(m)
    for (k,v) in bounds
        # Getting the bound value always as a t x r matrix, even if it has been expressed as a scalar or as a (temporal) vector
        v_matrix = Array{Float64}(undef,ntsteps,nreg)
        if (ndims(v[2]) == 0)
            v_matrix .= v[2]
        elseif (ndims(v[2]) == 1)
            [v_matrix[:,c] .= v[2] for c in 1:nreg]
        elseif (ndims(v[2]) == 2)
            v_matrix .= v[2]
        else
            error("Unknown constraint type $(k)")
        end 
        if (v[1] == "<=")
            scenboundseq = @constraint(m, [ti in tidx, ri in ridx], mvars[Symbol(k)][ti,ri] <= v_matrix[ti, ri])
        elseif (v[1] == ">=")
            scenboundseq = @constraint(m, [ti in tidx, ri in ridx], mvars[Symbol(k)][ti,ri] >= v_matrix[ti, ri])
        elseif (v[1] == "==")
            scenboundseq = @constraint(m, [ti in tidx, ri in ridx], mvars[Symbol(k)][ti,ri] == v_matrix[ti, ri])
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

    # Ex-post computation of global variables
    F_GHGABATE = [sum(value.(F_GHGABATE_R)[ti,ri] for ri in ridx) for ti in tidx]
    MIU        = [(1 - sum(value.(EIND_R)[ti,ri] for ri in ridx)/sum(sigma[ti,ri] * value.(YGROSS_R)[ti,ri] for ri in ridx)) for ti in tidx]
    DAMFRAC    = [ sum(value.(DAMFRAC_R)[ti,ri] for ri in ridx) / sum(value.(YGROSS_R)[ti,ri] for ri in ridx) for ti in tidx]

    # Other post-optimization computaitons
    rfactlong  =  collect([srf[ri] *(value.(CPC_R)[ti,ri]/value.(CPC_R[1,ri]))^(-elasmu[ri])*rr[ti,ri] for ti in tidx, ri in ridx])
    rlong      =  [-log(rfactlong[ti,ri]/srf[ri])/(tstep*ti) for ti in tidx,ri in ridx]
    rshort     =  [ ti == 1 ? missing : -log(rfactlong[ti,ri]/rfactlong[ti-1,ri])/tstep for ti in tidx, ri in ridx]

    scc        = collect(@. -1000*dual(eco2eq_r) /(0.00001+dual(cc_r)))
    ppm        = collect(value.(MAT) ./ 2.13)
    abaterat   = collect(@. value(ABATECOST)/value(Y))
    atfrac2020 = collect(@. (value(MAT)-mat0)/(value(CCATOT)+0.00001-cumemiss0))
    atfrac1765 = collect(@. (value(MAT)-mateq)/(value(CCATOT)+0.00001))
    forc_co2   = collect(@. fco22x*log(value(MAT)/mateq)/log(2))

    

    # Return results as named tuple
    return (solved=true, status=status, times=times, tidx=tidx, rfactlong=rfactlong,rlong=rlong, rshort=rshort, scc=scc, ppm=ppm, abaterat=abaterat, atfrac2020=atfrac2020, atfrac1765=atfrac1765, forc_co2=forc_co2,
      EIND_R=collect(value.(EIND_R)),EIND=collect(value.(EIND)),
      ECO2_R=collect(value.(ECO2_R)),ECO2=collect(value.(ECO2)),
      ECO2E_R=collect(value.(ECO2E_R)),ECO2E=collect(value.(ECO2E)),
      F_GHGABATE_R=collect(value.(F_GHGABATE_R)),F_GHGABATE=collect(value.(F_GHGABATE)),
      MIU_R=collect(value.(MIU_R)), MIU=collect(value.(MIU)),
      C=collect(value.(C)),C_R=collect(value.(C_R)),
      K=collect(value.(K)),K_R=collect(value.(K_R)),
      CPC_R=collect(value.(CPC_R)), CPC=collect(value.(CPC)),
      I_R=collect(value.(I_R)), I=collect(value.(I)),
      S_R=collect(value.(S_R)), S=collect(value.(S)),
      Y_R=collect(value.(Y_R)), Y=collect(value.(Y)),
      YGROSS_R=collect(value.(YGROSS_R)), YGROSS=collect(value.(YGROSS)),
      YNET_R=collect(value.(YNET_R)), YNET=collect(value.(YNET)),
      DAMAGES_R=collect(value.(DAMAGES_R)), DAMAGES=collect(value.(DAMAGES)),
      DAMFRAC_R=collect(value.(DAMFRAC_R)), DAMFRAC=collect(value.(DAMFRAC)),
      ABATECOST_R=collect(value.(ABATECOST_R)), ABATECOST=collect(value.(ABATECOST)),
      CCATOT=collect(value.(CCATOT)),
      PERIODU_R=collect(value.(PERIODU_R)), PERIODU=collect(value.(PERIODU)),
      CPRICE_R=value.(CPRICE_R),
      TOTPERIODU=collect(value.(TOTPERIODU)),
      UTILITY=value(UTILITY),
      FORC=collect(value.(FORC)),
      TATM=collect(value.(TATM)),
      TBOX1=collect(value.(TBOX1)), TBOX2=collect(value.(TBOX2)),
      RES0=collect(value.(RES0)), RES1=collect(value.(RES1)), RES3=collect(value.(RES3)),
      MAT=collect(value.(MAT)),
      CACC=collect(value.(CACC)),
      IRFT=collect(value.(IRFT)),
      ALPHA=collect(value.(ALPHA)),
      pars=pars, model=m)


end

function run_dice(;optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0, "max_iter" => 3000, ), bounds=Dict{String,Tuple{String,String}}(),kwargs...) 

    issubset(keys(kwargs), fieldnames(DICEParameters)) || error("Not all keywords are valid parameters.")
    pars = DICE2023(;kwargs...)   # Override the default parameter values with the keyword arguments
    ret = return run_dice(pars;optimizer=optimizer,bounds=bounds)
end