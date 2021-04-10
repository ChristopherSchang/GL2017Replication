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
plotly()
using Parameters
@unpack b_grid,c_pol,n_pol,b_pol,JD = gl


# From paper:   Figure 1
Y1 = 1. # !!! made up
plot(b_grid/(4*Y1), c_pol[2,:])
plot!(b_grid/(4*Y1), c_pol[8,:])
title!("consumption") 

plot(b_grid/(4*Y1), n_pol[2,:])
plot!(b_grid/(4*Y1), n_pol[8,:]) 
title!("labor supply") 

plot(b_grid/(4*Y1), b_pol[2,:])
plot!(b_grid/(4*Y1), b_pol[8,:]) 
title!("bond policy") 
 

