# DICEModel.jl

Documentation for `DICEModel.jl`, a Julia/JuMP port of the Nordhouse's DICE (Dynamic Integrated Climate-Economy model) model.


This package currently implements the DICE2023-b-4-3-10.gms gams version.

**This program and output is not the original Barrage/Nordhaus version, which is currently only [available in GAMS](https://bit.ly/3TwJ5nO).**


### Available documentation
- [`run_dice_scenario(scenario_name)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.run_dice_scenario-Tuple{String}): run one of the "official" 10 scenarios ([browse code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/Scenarios.jl));
- [`run_dice(;optimizer,bounds,kwargs...)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.run_dice-Tuple{}): run DICE with custom solver engine (and eventually options), custom variable constraints (bounds) or custom parameters ([browse code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/DICEModel.jl));
- [`Parameters`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.RawParameters): Available parameters to use with the `run_dice` function ([browse code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/DICEModel.jl))

### Results

In both cases the output (results) is a named tuple. Use `keys(results)` to find the available information (or just look at the source code) and `results.VARIABLEX` to obtain the values.

A summary of the main results, and a comparision with the official Barrage/Nordhaus version, is available [on this page](https://sylvaticus.github.io/DICEModel.jl/dev/results.html).

An Excel version of all the results is available here.

### Example

```julia
using Pkg
Pkg.activate(@__DIR__)
Pkg.add(["DICEModel","Plots"]) # run only once, then comment out
using DICEModel, Plots

# CB Optimal scenario...
res_cbopt    = run_dice_scenario("cbopt")

# Base scenario...
res_base    = run_dice_scenario("base")

# Paris "extended" scenario...
tidx = 1:81
# upper limit to emissions mitigation rate
miuup = @. min( 0.05 + 0.04*(tidx-1) - 0.01*max(0,tidx-5)  ,1.00) 
res_parisext = run_dice(miuup = miuup) # or simply: run_dice_scenario("parisext")

# Max 2 °C scenario...
res_t2c = run_dice(bounds = Dict("TATM"=>("<=",2.0))) # or simply: run_dice_scenario("t2c")

# Plots
times = res_cbopt.times

# CO2 emissions plot...
plot(times[1:11],res_cbopt.ECO2[1:11],ylim=(0,70), title="CO₂ emissions",ylabel="GtCO₂/yr",label="C/B optimal", color=:blue4, markershape=:circle, markercolor=:white)
plot!(times[1:11],res_base.ECO2[1:11], label="Base", colour=:goldenrod3, markershape=:circle, markercolor=:goldenrod3)
plot!(times[1:11],res_parisext.ECO2[1:11], label="Paris ext", colour=:red, linestyle=:dash)
plot!(times[1:11],res_t2c.ECO2[1:11], label="T < 2 °C", colour=:green, markershape=:cross, markercolor=:green)

# Carbon price plot...
plot(times[1:9],res_cbopt.CPRICE[1:9],ylim=(0,300), title="Carbon price",ylabel="2019\$ / t tCO₂",label="C/B optimal", color=:blue4, markershape=:circle, markercolor=:white)
plot!(times[1:9],res_base.CPRICE[1:9], label="Base", colour=:goldenrod3, markershape=:circle, markercolor=:goldenrod3)
plot!(times[1:9],res_parisext.CPRICE[1:9], label="Paris ext", colour=:red, linestyle=:dash)
plot!(times[1:9],res_t2c.CPRICE[1:9], label="T < 2 °C", colour=:green, markershape=:cross, markercolor=:green)
```

<img src="![assets/imgs/CO₂_emissions.png](https://github.com/sylvaticus/DICEModel.jl/blob/main/assets/imgs/CO%E2%82%82_emissions.png?raw=true)" width="400"/> <img src="![assets/imgs/Carbon_price.png](https://github.com/sylvaticus/DICEModel.jl/blob/main/assets/imgs/Carbon_price.png?raw=true)" width="400"/>



## Licence
The licence of the original GAMS code has never being specified. The Julia port itself (and only that) is MIT.



