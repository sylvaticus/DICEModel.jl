
#using Pkg 
#Pkg.activate(@__DIR__)
#Pkg.instantiate()

using Documenter
using DICEModel

push!(LOAD_PATH,"../src/")
makedocs(sitename="DICEModel.jl Documentation",
         pages = [
            "Index" => "index.md",
            "API" => "api.md",
            "Results" => "results.md"
         ],
         format = Documenter.HTML(prettyurls = false, assets=["assets/custom.js","assets/custom.css"], analytics = "G-0MDFVNSBYE"),
         warnonly = true,
         checkdocs=:none,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/sylvaticus/DICEModel.jl.git",
    devbranch = "main"
)
