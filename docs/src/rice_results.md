# Multi-regional (RICE) Results

This page reports preliminary results from running the package with multi-regional partitions of the world and compare them with the 1-region implementation (DICE).

We first run the model with _equal_ n regions (changing n) and then with the more realistic RICE-based partition (changing the weights).

**Warning** This is a first draft of the model results. Several anomalies are present, most likely due to errors in the code or in the parametrisation.


```@contents
Pages = ["rice_results.md"]
Depth = 4
```

We start by loading the required packages.

!!! details "Show code"
    ```@example mr
    using Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    #Pkg.add(["DICEModel", "CSV","DataFrames","Plots"]) # Run once and comment this back
    using DICEModel, Ipopt, JuMP, Plots, Test 
    nothing #hide
    ```

## Effect of multiple regions

We run `DICE2023_NREG` testing 1,2,4,6,8,10 and 12 regional partitions.

!!! details "Show code"
    ```@example mr
    nregs = [1,2,4,6,8,12]
    n_nregscen = length(nregs)
    results = Dict(["$(nregs[n])_regions" => run_dice(DICE2023_NREG(nregs[n]), optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5, "max_iter" => 4000, "tol"=> 5*10^-8, "acceptable_tol" =>5*10^-6)) for n in 1:n_nregscen])
    plot_attributes = Dict(
        "1_regions"    => (label="1 region (DICE)",   linestyle=:solid, colour=:red3, linewidth=:2),
        "2_regions"    => (label="2 regions",  linestyle=:solid, colour=:navy,linewidth=:3),
        "4_regions"    => (label="4 regions",  linestyle=:dash, colour=:midnightblue,linewidth=:3),
        "6_regions"    => (label="6 regions",  linestyle=:dot, colour=:indigo,linewidth=:3),
        "8_regions"    => (label="8 regions",  linestyle=:solid, colour=:blue3,linewidth=:1),
        "10_regions"   => (label="10 regions", linestyle=:dash, colour=:darkslategray,linewidth=:1),
        "12_regions"   => (label="12 regions", linestyle=:dot, colour=:royalblue2,linewidth=:1),
    )
    times  = results["1_regions"].times.+2020
    nothing #hide
    ```

### Emissions

!!! details "Show code"
    ```@example mr
    scenarios_plot=["$(nregs[n])_regions" for n in 1:n_nregscen]
    p = plot();
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:11] ,results[s].ECO2[1:11],ylim=(0,120), title="CO₂ emissions",ylabel="GtCO₂/yr";plot_attributes[s]...);
        else
            plot!(times[1:11] ,results[s].ECO2[1:11];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example mr
p #hide
```


### Temperatures

!!! details "Show code"
    ```@example mr
    p = plot();
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:17] ,results[s].TATM[1:17],ylim=(0.0,4.0), title="Global Temperatures",ylabel="°C from 1765";plot_attributes[s]...);
        else
            plot!(times[1:17] ,results[s].TATM[1:17];plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example mr
p #hide
```

### Emission Control rate

!!! details "Show code"
    ```@example mr
    for(i,s) in enumerate(scenarios_plot)
        if i == 1
            global p = plot(times[1:17] ,results[s].MIU[1:17] .* 100,ylim=(0.0,100), title="Emission control rate",ylabel="%";plot_attributes[s]...);
        else
            plot!(times[1:17] ,results[s].MIU[1:17] .* 100;plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example mr
p #hide
```

### Interregional variability

Testing that, with equal weighths, the regions behave the same:

```@example mr
@test results["2_regions"].ECO2_R[:,1] ≈ results["2_regions"].ECO2_R[:,2]
```

Testing a skewed welfare weights distribution, with more welfare weighth to the first region:

!!! details "Show code"
    ```@example mr
    results2r = run_dice(DICE2023_NREG(2,weights=[5,1]))
    p1 = plot(times[1:11] ,results2r.ECO2_R[1:11,1],ylim=(0,40), title="CO₂ emissions",ylabel="GtCO₂/yr", label="More utility weigthed region")
    plot!(p1,times[1:11] ,results2r.ECO2_R[1:11,2],ylim=(0,40), title="CO₂ emissions",ylabel="GtCO₂/yr", label="Less utility weigthed region")
    p2 = plot(times[1:17] ,results2r.MIU_R[1:17,1] .* 100, ylim=(0.0,100), title="Emission control rate",ylabel="%", label="More utility weigthed region");
    plot!(p2,times[1:17] ,results2r.MIU_R[1:17,2] .* 100, ylim=(0.0,100), title="Emission control rate",ylabel="%",label="Less utility weigthed region")
    nothing #hide
    ```

```@example mr
p1 #hide
```
```@example mr
p2 #hide
```


## RICE2023

Here we run `RICE2023`, a version of the multi-regional model with a parameterisation where regional data has a distribution from the RICE2020 model, but rescaled to total (world) data from DICE2023.

!!! details "Show code"
    ```@example mr
    rice_eq = run_dice(RICE2023(), optimizer=optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 5, "max_iter" => 3000, "tol"=> 5*10^-8, "acceptable_tol" =>5*10^-6))

    regions = rice_eq.pars.regions
    nreg = length(regions)
    nothing #hide
    ```

### Totals and comparision with DICE/DICE 12 region

!!! details "Show code"
    ```@example mr
    plot_attributes = Dict(
        "1_regions"    => (label="DICE",   linestyle=:solid, colour=:red3, linewidth=:2),
        "12_regions"  => (label="DICE 12 EQ REG",   linestyle=:dash, colour=:green3, linewidth=:2),
        "rice"  => (label="RICE 12 REG",   linestyle=:dot, colour=:blue3, linewidth=:2)
    )
    nothing #hide
    ```

#### Emissions

!!! details "Show code"
    ```@example mr
    p = plot(times[1:11] ,results["1_regions"].ECO2[1:11],ylim=(0,100), title="CO₂ emissions",ylabel="GtCO₂/yr";plot_attributes["1_regions"]...);
    plot!(times[1:11] ,results["12_regions"].ECO2[1:11];plot_attributes["12_regions"]...);
    plot!(times[1:11] ,rice_eq.ECO2[1:11];plot_attributes["rice"]...)
    nothing #hide
    ```
```@example mr
p #hide
```

#### Temperatures

!!! details "Show code"
    ```@example mr
    p = plot(times[1:11] ,results["1_regions"].TATM[1:11],ylim=(0.0,4.0), title="Global Temperatures",ylabel="°C from 1765";plot_attributes["1_regions"]...);
    plot!(times[1:11] ,results["12_regions"].TATM[1:11];plot_attributes["12_regions"]...);
    plot!(times[1:11] ,rice_eq.TATM[1:11];plot_attributes["rice"]...)
    nothing #hide
    ```
```@example mr
p #hide
```

#### Mean control rate

!!! details "Show code"
    ```@example mr
    p = plot(times[1:11] ,results["1_regions"].MIU[1:11] .* 100,ylim=(0.0,100), title="Emission control rate",ylabel="%";plot_attributes["1_regions"]...);
    plot!(times[1:11] ,results["12_regions"].MIU[1:11] .* 100;plot_attributes["12_regions"]...);
    plot!(times[1:11] ,rice_eq.MIU[1:11] .* 100;plot_attributes["rice"]...)
    nothing #hide
    ```
```@example mr
p #hide
```

### Regional Variance

#### Emissions

!!! details "Show code"
    ```@example mr
    plot_attributes = Dict(
        regions[1]    => (label=regions[1], linestyle=:solid, colour=:navy, linewidth=:2, fillalpha=0.8),
        regions[2]    => (label=regions[2], linestyle=:dash, colour=:navy, linewidth=:2, fillalpha=0.5),
        regions[3]    => (label=regions[3], linestyle=:dashdot, colour=:navy, linewidth=:2, fillalpha=0.3),
        regions[4]    => (label=regions[4], linestyle=:solid, colour=:chocolate2, linewidth=:2, fillalpha=0.8),
        regions[5]    => (label=regions[5], linestyle=:dash, colour=:chocolate2, linewidth=:2, fillalpha=0.5),
        regions[6]    => (label=regions[6], linestyle=:dashdot, colour=:chocolate2, linewidth=:2, fillalpha=0.3),
        regions[7]    => (label=regions[7], linestyle=:solid, colour=:red4, linewidth=:2, fillalpha=0.8),
        regions[8]    => (label=regions[8], linestyle=:dash, colour=:red4, linewidth=:2, fillalpha=0.5),
        regions[9]    => (label=regions[9], linestyle=:dashdot, colour=:red4, linewidth=:2, fillalpha=0.3),
        regions[10]    => (label=regions[10], linestyle=:solid, colour=:darkgreen, linewidth=:2, fillalpha=0.8),
        regions[11]    => (label=regions[11], linestyle=:dash, colour=:darkgreen, linewidth=:2, fillalpha=0.5),
        regions[12]    => (label=regions[12], linestyle=:dashdot, colour=:darkgreen, linewidth=:2, fillalpha=0.3),
    )

    plot_colours = hcat([plot_attributes[regions[i]].colour for i in 1:nreg]...)
    plot_alphas  = hcat([plot_attributes[regions[i]].fillalpha for i in 1:nreg]...)
    plot_rnames  = hcat([plot_attributes[regions[i]].label for i in 1:nreg]...)
    nothing #hide
    ```

### Emmissions

!!! details "Show code"
    ```@example mr
    p = areaplot(times[1:11], 100 .* rice_eq.ECO2_R[1:11,:] ./ sum(rice_eq.ECO2_R[1:11,:],dims=2), title= "Share of CO₂eq emissions", ylabel="%", seriescolor = plot_colours, fillalpha = plot_alphas, label = plot_rnames)
    nothing #hide
    ```
```@example mr
p #hide
```

### Control rates

!!! details "Show code"
    ```@example mr
    p = plot();
    for(i,s) in enumerate(regions)
        if i == 1
            global p = plot(times[1:11] ,rice_eq.MIU_R[1:11,i] .* 100,ylim=(0.0,100), title="Emission control rate",ylabel="%";plot_attributes[s]...);
        else
            plot!(times[1:11] ,rice_eq.MIU_R[1:11,i] .* 100;plot_attributes[s]...)
        end
    end
    nothing #hide
    ```
```@example mr
p #hide
```

### Welfare analysis (weights distribution)

Finally, we complement the default RICE213 output, that runs on equal welfare weights, with a couple of runs when we "favour" developed countries (giving them a higher weighth) or, at the opposite, less developed ones:

!!! details "Show code"
    ```@example mr
    w_rich  = [5,4,3,3,1,1,3,2,2,1,1.5,1]
    w_poor  = [1,1,1,1.5,2,2,3,3,2,5,5,5]

    rice_rich = run_dice(RICE2023(;weights=w_rich))
    rice_poor = run_dice(RICE2023(;weights=w_poor))

    rice_wanalysis = [rice_rich, rice_eq, rice_poor]
    ncols = 4
    nrows = 3

    plot_attributes = Dict(
        "rich"    => (label="Pref to developed countries",   linestyle=:solid, colour=:red3, linewidth=:2),
        "eq"  => (label="Equal preferences",   linestyle=:dash, colour=:green3, linewidth=:2),
        "poor"  => (label="Pref to underdeveloped countries",   linestyle=:dashdot, colour=:blue3, linewidth=:2)
    )
    nothing #hide
    ```

#### Emissions

##### Global

!!! details "Show code"
    ```@example mr
    p = plot(times[1:11] ,rice_rich.ECO2[1:11],ylim=(0,100), title="CO₂ emissions",ylabel="GtCO₂/yr";plot_attributes["rich"]...);
    plot!(times[1:11] ,rice_eq.ECO2[1:11];plot_attributes["eq"]...);
    plot!(times[1:11] ,rice_poor.ECO2[1:11];plot_attributes["poor"]...)
    nothing #hide
    ```
```@example mr
p #hide
```

##### Regional control rate

!!! details "Show code"
```@example mr
plots = Array{Plots.Plot,1}(undef,nreg)
for ri in 1:nreg
    c = ri%ncols == 0 ? ncols : ri%ncols
    r = Int(ceil(ri/ncols))
    ylabel = (c == 1) ? "%" : ""
    #legend = (ri == 1) ? :bottomright : nothing
    pr = plot(times[1:11],rice_wanalysis[1].MIU_R[1:11,ri] .* 100;ylabel=ylabel,ylim=(0.0,100.0),xlabel=regions[ri],legend=nothing,plot_attributes["rich"]...) 
    plot!(pr,times[1:11],rice_wanalysis[2].MIU_R[1:11,ri] .* 100;plot_attributes["eq"]...)
    plot!(pr,times[1:11],rice_wanalysis[3].MIU_R[1:11,ri] .* 100;plot_attributes["poor"]...)
    plots[ri] = pr
end
sizebottom= 3000
plot(plots...,layout=(nrows,ncols), size=(sizebottom/3,sizebottom/4))
nothing #hide
```
```@example mr
plot(plots...,layout=(nrows,ncols), size=(sizebottom/3,sizebottom/4)) #hide
```

#### Temperatures

!!! details "Show code"
    ```@example mr
    p = plot(times[1:11] ,rice_rich.TATM[1:11],ylim=(0.0,4.0), title="Global Temperatures",ylabel="°C from 1765";plot_attributes["rich"]...);
    plot!(times[1:11] ,rice_eq.TATM[1:11];plot_attributes["eq"]...);
    plot!(times[1:11] ,rice_poor.TATM[1:11];plot_attributes["poor"]...)
    nothing #hide
    ```
```@example mr
p #hide
```




