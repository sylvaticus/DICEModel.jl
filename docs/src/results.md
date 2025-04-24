# DICE Results

In this page we report the results from running the package with the official scenarios and compare the output of `DICEModel.jl` with the official GAMS output as given in the _DICE 2023 Introduction and User's Manual_ (v3.1.2, May 15, 2024)

The scenarios considered are:

- `cbopt`:    The C/B optimal scenario
- `t2c`:      The temperature constrained to max 2 °C scenario
- `t15c`:     The temperature constrained to max 1.5 °C scenario
- `altdam`:   The alternative damage scenario
- `parisext`: The Paris extended scenario
- `base`:     The base (current policies) scenario
- `r5`:       The scenario with discount rate of 5%
- `r4`:       The scenario with discount rate of 4%
- `r3`:       The scenario with discount rate of 3%
- `r2`:       The scenario with discount rate of 2%
- `r1`:       The scenario with discount rate of 1%

```@contents
Pages = ["results.md"]
Depth = 4
```

We start by loading the required packages.

!!! details "Show code"
    ```@example r
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    #Pkg.add(["DICEModel", "CSV","DataFrames","Plots"]) # Run once and comment this back
    using DICEModel, CSV, DataFrames, Plots, XLSX
    nothing #hide
    ```

We can now define the scenarios to run and their specific graphical attributed when plotted. 
Finally we call the function `run_dice_scenario` with each of them and save the output in the `results` dictionary (keyed by the scenario name):

!!! details "Show code"
    ```@example r
    # Scenarios to run:
    scenarios = ["cbopt","t2c","t15c","altdam","parisext","base","r5","r4","r3","r2","r1"]
    plot_attributes = Dict(
        "cbopt"    => (label="C/B optimal", linestyle=:solid, colour=:red, markershape=:circle, markercolor=:white, markerstrokecolor=:red),  
        "t2c"      => (label="T < 2 °C", linestyle=:dash, colour=:green, markershape=:cross, markercolor=:green, markerstrokecolor=:green),     
        "t15c"     => (label="T < 1.5 °C", linestyle=:dot, colour=:lime, markershape=:dtriangle, markercolor=:lime, markerstrokecolor=:lime),     
        "altdam"   => (label="Alt. damage", linestyle=:solid, colour=:darkgoldenrod4, markershape=:none, markercolor=:darkgoldenrod4, markerstrokecolor=:darkgoldenrod4),   
        "parisext" => (label="Paris ext", linestyle=:dash, colour=:darkorange, markershape=:cross, markercolor=:darkorange, markerstrokecolor=:darkorange), 
        "base"     => (label="Base", linestyle=:solid, colour=:yellow, markershape=:utriangle, markercolor=:yellow, markerstrokecolor=:yellow),     
        "r5"       => (label="R = 5%", linestyle=:solid, colour=:darkblue, markershape=:diamond, markercolor=:darkblue, markerstrokecolor=:darkblue),       
        "r4"       => (label="R = 4%", linestyle=:dash, colour=:blue, markershape=:diamond, markercolor=:blue, markerstrokecolor=:blue),       
        "r3"       => (label="R = 3%", linestyle=:solid, colour=:deepskyblue4, markershape=:diamond, markercolor=:deepskyblue4, markerstrokecolor=:deepskyblue4),     
        "r2"       => (label="R = 2%", linestyle=:dash, colour=:steelblue, markershape=:diamond, markercolor=:steelblue, markerstrokecolor=:steelblue),        
        "r1"       => (label="R = 1%", linestyle=:solid, colour=:lightskyblue, markershape=:diamond, markercolor=:lightskyblue, markerstrokecolor=:lightskyblue)      
    )
    # Run the scenarios and collect their results in the "results" dictionary
    results = Dict([s => run_dice_scenario(s) for s in scenarios])
    times  = results["cbopt"].times.+2020
    nothing #hide
    ```

## Main results and comparison with the official GAMS version

For the variables `ECO2`, `ppm`, `TATM`, `MIU`, `CPRICE` and `scc`, we show the official values, the ones computed with DICEModel.jl and we plot these values.

---

### `ECO2` : CO₂ emissions [GtCO₂/yr]

#### Official GAMS outputs:

!!! details "Show code"
    ```@example r
    # Data from the DICE2023 manual
    ECO2_gams = CSV.read(IOBuffer("""
    scen      2020  2025  2050  2100     
    cbopt     42.9  42.9  37.1  15.9 
    t2c       42.9  42.9  27.2   1.2
    t15c      42.9  13.1   5.7   0.0
    altdam    42.9  42.7  20.9   0.0  
    parisext  42.9  43.3  44.4  42.3    
    base      42.9  44.9  54.6  75.7 
    r5        42.8  42.5  42.2  37.6
    r4        44.1  43.9  39.3  28.9
    r3        45.6  45.3  33.5  15.4
    r2        46.8  46.7  22.5   0.0
    r1        46.8  46.9  19.2   0.0
    """), DataFrame, delim=" ", ignorerepeated=true)
    nothing #hide
    ```

```@example r
ECO2_gams  #hide
```

#### DICEModel.jl output

!!! details "Show code"
    ```@example r
    ECO2 = DataFrame([
        scenarios,
        [results[s].ECO2[1] for s in scenarios],
        [results[s].ECO2[2] for s in scenarios],
        [results[s].ECO2[7] for s in scenarios],
        [results[s].ECO2[17] for s in scenarios],
        ],
        ["scen","2020","2025","2050","2100"]
    )
    nothing #hide
    ```

```@example r
ECO2 #hide
```

#### Differences:

!!! details "Show code"
    ```@example r
    ECO2_diff = copy(ECO2_gams);
    ECO2_diff[:,2:end] .= ECO2_gams[:,2:end] .- ECO2[:,2:end]
    nothing #hide
    ```

```@example r
ECO2_diff #hide
```

There are only two minor differences in `r5` for 2100 and, above all, for `r1` in 2020.

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["base","cbopt","parisext","t2c","r2","r1"]
    p = plot()
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:11] ,results[s].ECO2[1:11],ylim=(0,80), title="CO₂ emissions",ylabel="GtCO₂/yr";plot_attributes[s]...);
        else
            plot!(times[1:11] ,results[s].ECO2[1:11];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

### `ppm`: CO₂ concentration [ppm]

#### Official GAMS outputs:

!!! details "Show code"
    ```@example r
    # Data from the DICE2023 manual
    ppm_gams = CSV.read(IOBuffer("""
    scen       2020  2025  2050  2100   2150
    cbopt     416.2 429.9 487.8 569.2  497.9 
    t2c       416.2 429.9 474.7 474.7  437.9      
    altdam    416.2 429.8 466.7 458.5  401.0  
    parisext  416.2 430.1 501.3 652.5  763.5    
    base      416.2 430.9 517.7 774.9 1144.0 
    r5        416.2 429.7 495.7 635.5  671.6
    r4        416.2 430.4 491.8 605.3  592.0
    r3        416.2 431.2 484.0 555.1  494.2
    r2        416.2 431.9 473.8 484.9  419.3
    r1        416.2 432.0 473.3 449.5  389.3
    """), DataFrame, delim=" ", ignorerepeated=true)
    nothing #hide
    ```
```@example r
ppm_gams #hide
```

#### DICEModel.jl output

!!! details "Show code"
    ```@example r
    scenarios_ppm = ["cbopt","t2c","altdam","parisext","base","r5","r4","r3","r2","r1"]
    ppm = DataFrame([
        scenarios_ppm,
        [results[s].ppm[1] for s in scenarios_ppm],
        [results[s].ppm[2] for s in scenarios_ppm],
        [results[s].ppm[7] for s in scenarios_ppm],
        [results[s].ppm[17] for s in scenarios_ppm],
        [results[s].ppm[27] for s in scenarios_ppm],
        ],
        ["scen","2020","2025","2050","2100","2150"]
    )
    nothing #hide
    ```

```@example r
ppm #hide
```

#### Differences:

!!! details "Show code"
    ```@example r
    ppm_diff = copy(ppm_gams);
    ppm_diff[:,2:end] .= ppm_gams[:,2:end] .- ppm[:,2:end]
    nothing #hide
    ```

```@example r
ppm_diff #hide
```

Again, the only minor difference is for `r5` in 2150.

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["base","cbopt","parisext","t2c","r2","r1"]
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:17] ,results[s].ppm[1:17],ylim=(400,800), title="CO₂ concentrations",ylabel="ppm";plot_attributes[s]...);
        else
            plot!(times[1:17] ,results[s].ppm[1:17];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

### `TATM`: Global Temperatures [°C diff since 1765]

#### Official GAMS outputs:

!!! details "Show code"
    ```@example r
    # Data from the DICE2023 manual
    TATM_gams = CSV.read(IOBuffer("""
    scen       2020  2025  2050  2100  2150
    cbopt      1.25  1.42  1.92  2.58  2.29
    t2c        1.25  1.42  1.85  2.00  1.86     
    altdam     1.25  1.42  1.81  1.89  1.58 
    parisext   1.25  1.43  2.01  3.00  3.61   
    base       1.25  1.43  2.10  3.55  4.91
    r5         1.25  1.42  1.97  2.93  3.24
    r4         1.25  1.43  1.95  2.77  2.84
    r3         1.25  1.43  1.90  2.49  2.26
    r2         1.25  1.43  1.84  2.07  1.73
    r1         1.25  1.43  1.84  1.81  1.49
    """), DataFrame, delim=" ", ignorerepeated=true)
    nothing #hide
    ```
```@example r
TATM_gams #hide
```

#### DICEModel.jl output

!!! details "Show code"
    ```@example r
    scenarios_TATM = ["cbopt","t2c","altdam","parisext","base","r5","r4","r3","r2","r1"]
    TATM = DataFrame([
        scenarios_TATM,
        [results[s].TATM[1] for s in scenarios_TATM],
        [results[s].TATM[2] for s in scenarios_TATM],
        [results[s].TATM[7] for s in scenarios_TATM],
        [results[s].TATM[17] for s in scenarios_TATM],
        [results[s].TATM[27] for s in scenarios_TATM],
        ],
        ["scen","2020","2025","2050","2100","2150"]
    )
    nothing #hide
    ```
```@example r
TATM #hide
```

#### Differences:

!!! details "Show code"
    ```@example r
    TATM_diff = copy(TATM_gams);
    TATM_diff[:,2:end] .= TATM_gams[:,2:end] .- TATM[:,2:end]
    nothing #hide
    ```
```@example r
TATM_diff #hide
```
Again, the only minor difference is for `r5` in 2150 (3.25 instead of 3.24)

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["base","cbopt","parisext","t2c","r4","r2"]
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:17] ,results[s].TATM[1:17],ylim=(0.0,4.0), title="Global Temperatures",ylabel="°C from 1765";plot_attributes[s]...);
        else
            plot!(times[1:17] ,results[s].TATM[1:17];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

### `MIU`: Emission control rate [%]

#### Official GAMS outputs:

!!! details "Show code"
    ```@example r
    # Data from the DICE2023 manual
    MIU_gams = CSV.read(IOBuffer("""
    scen    2020 2030 2040 2050 2060 2100
    cbopt      5   24   31   39   46   84
    t2c        5   24   42   55   69   99   
    altdam     5   24   48   65   76  100
    parisext   5   13   21   27   33   57 
    base       5    6    8   10   12   22
    r5         5   19   23   29   34   60
    r4         5   24   29   36   42   70
    r3         5   24   39   47   54   85
    r2         5   24   48   66   73  100
    r1         5   24   48   72   90  100
    """), DataFrame, delim=" ", ignorerepeated=true)
    nothing #hide
    ```
```@example r
MIU_gams #hide
```

#### DICEModel.jl output

!!! details "Show code"
    ```@example r
    scenarios_MIU = ["cbopt","t2c","altdam","parisext","base","r5","r4","r3","r2","r1"]
    MIU = DataFrame([
        scenarios_MIU,
        [results[s].MIU[1] for s in scenarios_MIU] .* 100,
        [results[s].MIU[3] for s in scenarios_MIU] .* 100,
        [results[s].MIU[5] for s in scenarios_MIU] .* 100,
        [results[s].MIU[7] for s in scenarios_MIU] .* 100,
        [results[s].MIU[9] for s in scenarios_MIU] .* 100,
        [results[s].MIU[17] for s in scenarios_MIU] .* 100,
        ],
        ["scen","2020","2030","2040","2050","2060","2100"]
    )
    nothing #hide
    ```
```@example r
MIU #hide
```

#### Differences:

!!! details "Show code"
    ```@example r
    MIU_diff = copy(MIU);
    MIU_diff[:,2:end] .= MIU_gams[:,2:end] .- MIU[:,2:end]
    nothing #hide
    ```
```@example r
MIU_diff #hide
```

No differences here !

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["cbopt","parisext","t2c","altdam"]
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:17] ,results[s].MIU[1:17] .* 100,ylim=(0.0,100), title="Emission control rate",ylabel="%";plot_attributes[s]...);
        else
            plot!(times[1:17] ,results[s].MIU[1:17] .* 100;plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

### `CPRICE`: Carbon price [2019$ / tCO₂]

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["cbopt","parisext","t2c","r3"]
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:9] ,results[s].CPRICE[1:9],ylim=(0,350), title="Carbon price",ylabel="2019\$ / tCO₂";plot_attributes[s]...);
        else
            plot!(times[1:9] ,results[s].CPRICE[1:9];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

### `scc`: Social cost of carbon [2019\$ / tCO₂]

#### Official GAMS outputs:

!!! details "Show code"
    ```@example r
    # Data from the DICE2023 manual
    scc_gams = CSV.read(IOBuffer("""
    scen       2020  2025  2050
    cbopt        50    59   125
    t2c          75    89   213
    t15c       3557  4185 16552          
    altdam      124   146   281
    parisext     61    72   159     
    base         66    78   175
    r5           32    37    74
    r4           49    58   107
    r3           87   102   172         
    r2          176   207   302
    r1          485   571   695
    """), DataFrame, delim=" ", ignorerepeated=true)
    nothing #hide
    ```
```@example r
scc_gams #hide
```

#### DICEModel.jl output

!!! details "Show code"
    ```@example r
    scenarios_scc = ["cbopt","t2c","t15c", "altdam","parisext","base","r5","r4","r3","r2","r1"]
    scc = DataFrame([
        scenarios_scc,
        [results[s].scc[1] for s in scenarios_scc],
        [results[s].scc[2] for s in scenarios_scc],
        [results[s].scc[7] for s in scenarios_scc],
        ],
        ["scen","2020","2025","2050"]
    )
    nothing #hide
    ```
```@example r
scc #hide
```

#### Differences:

!!! details "Show code"
    ```@example r
    scc_diff = copy(scc);
    scc_diff[:,2:end] .= scc_gams[:,2:end] .- scc[:,2:end]
    nothing #hide
    ```
```@example r
scc_diff #hide
```

This is the only part that still needs to be checked, as there are significant differences for the base year (2020). For the other years `DICEModel.jl` provides identical results (up to the approximation of the data published) than the official GAMS version.

#### Plot:

!!! details "Show code"
    ```@example r
    scenarios_plot=["base","cbopt","t2c","altdam","r4","r2"]
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:7] ,results[s].scc[1:7],ylim=(0,400), title="Social cost of carbon",ylabel="2019\$ / tCO₂";plot_attributes[s]...);
        else
            plot!(times[1:7] ,results[s].scc[1:7];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example r
p #hide
```

---

## Detailed model output

The whole model output (variables and post-processing computed values) can be retrieved in the Excel file below:

!!! details "Show code"
    ```@example r
    keys(results["cbopt"]) 
    out_vars = [v for v in keys(results["cbopt"]) if typeof(results["cbopt"][v]) <: Vector{<: Number}]

    all_results = DataFrame(
        scenario = String[],
        variable = String[],
        year     = Int64[],
        value    = Float64[],
    )

    for s in scenarios, v in out_vars, ti in 1:length(times)
        push!(all_results,[s,string(v),times[ti],results[s][v][ti]])
    end
    all_results = unstack(all_results,"year","value")

    XLSX.writetable("DICEModelDetailedResults.xlsx",
    all_results=(collect(DataFrames.eachcol(all_results)),DataFrames.names(all_results))
    )
    nothing #hide
    ```

[DICEModelDetailedResults.xlsx](DICEModelDetailedResults.xlsx)
