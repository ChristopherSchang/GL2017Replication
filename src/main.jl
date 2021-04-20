using Revise
include("GL2017Replication.jl")
using .GL2017Replication

# -----------------------------------------------------------------------------

# Load model structure
gl = ModelGL() 

# solve model with given parameters
initilize!(gl)
EGM!(gl)
gl.JD = compute_distribution!(gl)
aggregate!(gl)

# calibrate to target values 
calibrate!(gl)

 
# -----------------------------------------------------------------------------
# plot results

using Plots
plotly(linewidth = 2.)
using Parameters
@unpack b_grid,c_pol,n_pol,b_pol,JD, ϕ,pr = gl

# From paper:   Figure I

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
  

# From paper:   Figure IV
using Statistics 
bond_distribution = dropdims( pr'*JD[:,b_grid .>= -ϕ] ,dims=1)  
plot(b_grid[b_grid .>= -ϕ]/(4*Y1),  bond_distribution    )
title!("bond distribution") 
xaxis!([-ϕ,12.5])
 