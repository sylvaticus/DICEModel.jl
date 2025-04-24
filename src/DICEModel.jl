"""
   DICEModel

Implementation of the DICE/RICE models

# Notes:
- Optimization based on DICE2023-b-4-3-10.gms and included files (Nonco2-b-4-3-1.gms and FAIR-beta-4-3-1.gms)
- Variable casing has been harmonized that all parameters and post-optimization computation have lower cases, and all optimization variables have upper case.
- DICE has been generalized to work with multiple regions, and it is present in particular a version that uses the structure of DICE 2023 and the regional distribution of initial values (production, emissions, capital, population...) of RICE2020.
"""
module DICEModel

export run_dice, run_dice_scenario, DICEParameters, DICE2023, DICE2023_NREG, RICE2023
import Random
#using PrecompileTools
using DocStringExtensions # just for precompilation and documentation 
using JuMP, Ipopt

include("Utilities.jl")       # Utility functions such `fields_to_vars` and `scaleweights`
include("Parameters.jl")      # Contain the definition of the parameters struct and several functions with predefined "defaults" (DICE2023, RICE2023...)
include("CoreModel.jl")       # Contains `run_dice`, the low level function that runs the optimizations 
include("Scenarios.jl")       # Implementation of `run_dice_scenario` with the "official" DICE2023 named scenarios
include("Precompilation.jl")  # Precompilation stuff needed for performances



end # module DICEModel