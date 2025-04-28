# DICEModel.jl

Documentation for [`DICEModel.jl`](https://github.com/sylvaticus/DICEModel.jl), a Julia/JuMP port of the Nordhaus' DICE (Dynamic Integrated Climate-Economy model) and Nordhaus/Yang RICE (Regional Integrated Climate-Economy model) models.

This package currently implements "exactly" (in structure, data and hence output) the DICE2023-b-4-3-10.gms gams version and generalize it introducing a regional dimension that can be used to implement RICE-like models.

**This program and output is not the original Barrage/Nordhaus version, which is currently only [available in GAMS](https://bit.ly/3TwJ5nO).**

While `DICEModel.jl` is _implemented_ in Julia it can easily be used in Python or R using the [`JuliaCall` (Python)](https://github.com/JuliaPy/PythonCall.jl) and [`JuliaCall` (R)](https://cran.r-project.org/web/packages/JuliaCall/index.html) packages respectively. 


### Available documentation
- [`run_dice_scenario(scenario_name)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.run_dice_scenario-Tuple{String}): run one of the "official" 10 scenarios ([browse code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/Scenarios.jl));
- [`run_dice(;optimizer,bounds,kwargs...)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.run_dice-Tuple{}): run DICE with custom solver engine (and eventually options), custom variable constraints (bounds) or custom parameters ([browse code](https://github.com/sylvaticus/DICEModel.jl/blob/main/src/CoreModel.jl));
- [`DICEParameters`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.DICEParameters): Available parameters to use with the `run_dice` function
- [`DICE2023(;kwargs...)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.DICE2023-Tuple{}): Instantiate a `DICEParameters` struct with defaults parameters to DICE2023 (single region)
- [`DICE2023_NREG(n;kwargs...)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.DICE2023_NREG): Build parameters for a DICE2023 world partitioned in _n_ equal regions (unless parameters are overrided)
- [`RICE2023(;kwargs...)`](https://sylvaticus.github.io/DICEModel.jl/dev/api.html#DICEModel.RICE2023-Tuple{}): Build parameters calibrated to have the DICE2023 world totals and the 12-regions RICE2020 regional distribution.

### Results

In both cases the output (results) is a named tuple. Use `keys(results)` to find the available information (or just look at the source code) and `results.VARIABLEX` to obtain the values.

A summary of the main results, and a comparision with the official Barrage/Nordhaus DICE2023 version, is available [on this page](https://sylvaticus.github.io/DICEModel.jl/dev/results.html).

An Excel version of all the results for DICE2023 is available [here](DICEModelDetailedResults.xlsx).

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

```@raw html
<img src="https://github.com/sylvaticus/DICEModel.jl/blob/main/assets/imgs/CO%E2%82%82_emissions.png?raw=true" width="300"/><img src="https://github.com/sylvaticus/DICEModel.jl/blob/main/assets/imgs/Carbon_price.png?raw=true" width="300"/>
```

## Licence
The licence of the original GAMS code has never being specified. The Julia port itself (and only that) is MIT.



