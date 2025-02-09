
#using Pkg 
#Pkg.activate(@__DIR__)
#Pkg.instantiate()

using Documenter
using JuliaDICE

push!(LOAD_PATH,"../src/")
makedocs(sitename="JuliaDICE.jl Documentation",
         pages = [
            "Index" => "index.md",
            "API" => "api.md",
         ],
         format = Documenter.HTML(prettyurls = false),
         warnonly = true,
         checkdocs=:none,

)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/sylvaticus/JuliaDice.jl.git",
    devbranch = "main"
)
