using Revise
include("GL2017Replication.jl")
using .GL2017Replication

# Load model structure
gl = ModelGL() 

# solve model
EGM!(gl)

# Compute invariant joint distribution
compute_distribution!(gl)


# plot results
using Plots
plotly(linewidth = 2.)
using Parameters
@unpack b_grid,c_pol,n_pol,b_pol,JD, ϕ₁ = gl
ϕ = ϕ₁

# From paper:   Figure 1
Y1 = 0.4174 # !!! made up
plot(b_grid[b_grid .>= -ϕ]/(4*Y1), c_pol[2,b_grid .>= -ϕ],label = "θ = 2")
plot!(b_grid[b_grid .>= -ϕ]/(4*Y1), c_pol[8,b_grid .>= -ϕ],
        linestyle = :dash,label = "θ = 8")
title!("consumption") 
xaxis!([-ϕ,12.5])

plot(b_grid[b_grid .>= -ϕ]/(4*Y1), n_pol[2,b_grid .>= -ϕ],label =  "θ = 2")
plot!(b_grid[b_grid .>= -ϕ]/(4*Y1), n_pol[8,b_grid .>= -ϕ],
        linestyle = :dash,label = "θ = 8") 
title!("labor supply") 
xaxis!([-ϕ,12.5])

plot(b_grid[b_grid .>= -ϕ]/(4*Y1), b_pol[2,b_grid .>= -ϕ])
plot!(b_grid[b_grid .>= -ϕ]/(4*Y1), b_pol[8,b_grid .>= -ϕ]) 
title!("bond policy") 
xaxis!([-ϕ,12.5])


