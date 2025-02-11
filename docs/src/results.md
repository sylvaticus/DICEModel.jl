# Main results

We report here some results from running the codepackage with the official scenarios and compare the output of `DICEModel.jl` with the official GAMS output as given in the _DICE 2023 Introduction and User's Manual_ (v3.1.2, May 15, 2024)

The scenarios we run here are:

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

We start by loading the required packages.

```@example r
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
#Pkg.add(["DICEModel", "CSV","DataFrames","Plots"]) # Run once and comment this back
using DICEModel, CSV, DataFrames, Plots
```


!!! details "Show code"
    ```@example 1
    a = [1,2,3]
    b = vcat(a,1)
    ```




```@example r
# Scenarios to run
scenarios = ["cbopt","t2c","t15c","altdam","parisext","base","r5","r4","r3","r2","r1"]
# Run them and collect the results in the "results" dictionary
results = Dict([s => run_dice_scenario(s) for s in scenarios])
times  = results["cbopt"].times.+2020
```

cbopt
t2c
t15c
altdam
parisext
base
r5
r4
r3
r2
r1

### ECO2 : CO₂ emissions [GtCO₂/yr]

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
plot(times[1:11] .+ 2020,results["cbopt"].ECO2[1:11],ylim=(0,80), title="CO₂ emissions",ylabel="GtCO₂/yr",label="C/B optimal", color=:blue4, markershape=:circle, markercolor=:white);
plot!(times[1:11] .+ 2020,results["t2c"].ECO2[1:11], label="T < 2 °C", colour=:green, markershape=:cross, markercolor=:green);
plot!(times[1:11] .+ 2020,results["parisext"].ECO2[1:11], label="Paris ext", colour=:red, linestyle=:dash)
plot!(times[1:11] .+ 2020,results["base"].ECO2[1:11], label="Base", colour=:goldenrod3, markershape=:circle, markercolor=:goldenrod3)
```