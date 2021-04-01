module GL2017Replication
using BenchmarkTools


export hello, domath, solve_steady_state!

# Import scripts
include("model_structures.jl")



""" Function to solve steady-state objects.  
Inputs
    p::ModelParameters          :   Structure holding parameter Values
    grid::ModelSolutions        :   Structure for grids
    solution::ModelSolutions    :   Structure for solutions
"""
function solve_steady_state!(p::ModelParameters,grid::ModelGrid,solution::ModelSolutions)
    #do stuff

end


"""
    hello(who::String)

Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)

Return `x + 5`.
"""
domath(x::Number) = x + 5

end
