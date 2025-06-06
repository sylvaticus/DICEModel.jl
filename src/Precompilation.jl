"""
Part of [DICEModel](https://github.com/sylvaticus/DICEModel.jl). Licence is MIT.

This file runs the default optimization in order to precompile the code. It is automatically executed when the package is installed (built).
"""

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to the DICEModel package or not (on Julia 1.8 and higher)
        run_dice()
    end
end