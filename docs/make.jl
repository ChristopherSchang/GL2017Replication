push!(LOAD_PATH,"../src/")
using Documenter, GL2017Replication

makedocs(modules = [GL2017Replication], sitename = "GL2017Replication.jl"
        ,pages = [
                "Index" => "index.md"]
                )

deploydocs(repo = "github.com/ChristopherSchang/GL2017Replication.jl.git", devbranch = "main")
