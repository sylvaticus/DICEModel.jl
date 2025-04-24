"""
Part of [DICEModel](https://github.com/sylvaticus/DICEModel.jl). Licence is MIT.

This file contains some utility functions.
"""



"""
    @fields_to_vars(t,x)

Utility macro to convert struct fields to local variables (for readibility, so that we can write `parameterx` instead of using everywhere `p.parameterx`).
"""
macro fields_to_vars(t::Symbol, x)
    type = Core.eval(__module__, t)
        if !isstructtype(type)
            throw(ArgumentError("@fieldvars only takes struct types, not $type."))
        end 
    esc(:( (; $(fieldnames(type)...)) = $x::$type ))
end

"""
    scaleweights(w)

Scale a vector of weights such that their sum is 1
"""
scaleweights(w) = w ./ sum(w)

