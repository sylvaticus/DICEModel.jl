# JuliaDICE
Port in Julia/JuMP of the Nordhouse's DICE (Dynamic Integrated Climate-Economy model) model. Currently v2023


The model currently implements the DICE2023-b-4-3-10.gms gams version.


**"This program and output is not the original Barrage/Nordhaus version, which is currently only available [in GAMS](https://bit.ly/3TwJ5nO)."**

The licence of the original GAMS code in unknown. The Julia port is open source, as are all the library used.

Two functions are provided: `run_dice(scnario_name)` to run one of the official 10 scenarios, or `run_dice(;optimizer,bounds,kwargs...)` to specific custom solver engine (and eventually its option), equality or inequality constraints on any variable or to override some parameters. This function is the one internally used by the `run_dice(scnario_name)` function. 