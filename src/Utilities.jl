""" Utility functions for the DICEModel.jl package"""

"""
    @fields_to_vars(t,x)

Utility macro to convert struct fields to local variables (for readibility, so that we can write `parx` instead of using everywhere `p.parx`).
"""
macro fields_to_vars(t::Symbol, x)
    type = Core.eval(__module__, t)
        if !isstructtype(type)
            throw(ArgumentError("@fieldvars only takes struct types, not $type."))
        end 
    esc(:( (; $(fieldnames(type)...)) = $x::$type ))
end

scaleweights(w) = w ./ sum(w)

# """ LogSumExp for efficiently computing log(sum(exp.(x))) """
# lse(x) = maximum(x)+log(sum(exp.(x .- maximum(x))))
# """softmax (x; β=1) \n\n The input x is a vector. Return a PMF"""
# softmax(x; β=one.(x)) = exp.((β .* x) .- lse(β .* x))  # efficient implementation of softmax(x)  = exp.(x) ./  sum(exp.(x))
# softmax(x, β) = softmax(x, β=β)