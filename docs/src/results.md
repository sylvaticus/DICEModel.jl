# Main results

We report here some results from running the codepackage with the official scenarios and compare the output of `DICEModel.jl` with the official GAMS output as given in the _DICE 2023 Introduction and User's Manual_ (v3.1.2, May 15, 2024)

The scenarios we run here are:

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

We start by loading the required packages.

```@example r
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
#Pkg.add(["DICEModel", "CSV","DataFrames","Plots"]) # Run once and comment this back
using DICEModel, CSV, DataFrames, Plots
```

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

# Run them and collect the results in the "results" dictionary
results = Dict([s => run_dice_scenario(s) for s in scenarios])
times  = results["cbopt"].times.+2020

base_attributes = (label="Base", markershape=:utriangle, markercolor=:red, colour=:red)
```

---

### `ECO2` : CO₂ emissions [GtCO₂/yr]

#### Official GAMS outputs:

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
```

#### DICEModel.jl output

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
```
#### Differences:

```@example r
ECO2_diff = copy(ECO2_gams);
ECO2_diff[:,2:end] .= ECO2_gams[:,2:end] .- ECO2[:,2:end]
ECO2_diff
```

As you can see, there are only two minor differences in `r5` for 2100 and, above all, for `r1` in 2020.


#### Plot:

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
p
```

---

### `ppm`: CO₂ concentration [ppm]

#### Official GAMS outputs:

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
```

#### DICEModel.jl output

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
```
#### Differences:

```@example r
ppm_diff = copy(ppm_gams);
ppm_diff[:,2:end] .= ppm_gams[:,2:end] .- ppm[:,2:end]
ppm_diff
```

As you can see, there are only two minor differences in `r5` for 2100 and, above all, for `r1` in 2020.

#### Plot:

```@example r
scenarios_plot=["base","cbopt","parisext","t2c","r2","r1"]
for(i,s) in enumerate(scenarios_plot)
    if i == 1
        global p = plot(times[1:17] ,results[s].ppm[1:17],ylim=(400,800), title="CO₂ concentrations",ylabel="ppm";plot_attributes[s]...);
    else
        plot!(times[1:17] ,results[s].ppm[1:17];plot_attributes[s]...)
    end
end
p
```

