
using Parameters
 
@with_kw mutable struct ModelParameters

    # Deep parameters
    β::Float64               = 0.96  # discount rate
 
end

@with_kw mutable struct ModelGrid
    
    # grids
    income::Vector{Float64, 1}           = ones(1) 
 
end
 
@with_kw mutable struct ModelSolutions
    
    # Value function
    value::Array{Float64, 3}             = ones(1) 
 
end

 
 