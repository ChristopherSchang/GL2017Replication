push!(LOAD_PATH,"../src/")
using Documenter, GL2017Replication

makedocs(modules = [GL2017Replication], sitename = "GL2017Replication.jl"
        ,pages = [
                "Home" => "index.md",
                "Steady State" => "steadystate.md",
                "Transition" => "transition.md",
                "Output" => "output.md",
                "Functions" => "functions.md",
                ]
                )

deploydocs(repo = "github.com/ChristopherSchang/GL2017Replication.jl.git", devbranch = "main")
