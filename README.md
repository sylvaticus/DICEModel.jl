# JuliaDICE
Port in Julia/JuMP of the Nordhouse's DICE (Dynamic Integrated Climate-Economy model) model. Currently v2023


This package currently implements the DICE2023-b-4-3-10.gms gams version.

[![Build status (Github Actions)](https://github.com/sylvaticus/JuliaDICE.jl/workflows/CI/badge.svg)](https://github.com/sylvaticus/JuliaDICE.jl/actions)
[![codecov.io](http://codecov.io/github/sylvaticus/JuliaDICE.jl/coverage.svg?branch=main)](http://codecov.io/github/sylvaticus/JuliaDICE.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sylvaticus.github.io/JuliaDICE.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sylvaticus.github.io/JuliaDICE.jl/dev)


**"This program and output is not the original Barrage/Nordhaus version, which is currently only available [in GAMS](https://bit.ly/3TwJ5nO)."**

The licence of the original GAMS code has never being specified. The Julia port itself (and only that) is MIT.

Two functions are provided: `run_dice(scnario_name)` to run one of the official 10 scenarios, or `run_dice(;optimizer,bounds,kwargs...)` to specific custom solver engine (and eventually its option), equality or inequality constraints on any variable or to override some parameters. This function is the one internally used by the `run_dice(scnario_name)` function. 