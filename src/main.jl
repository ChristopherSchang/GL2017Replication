using Revise, JLD
include("GL2017Replication.jl")
using .GL2017Replication

# -----------------------------------------------------------------------------

# 1) Load model structure
gl = ModelGL()

# 2) solve model with default parameters
initilize!(gl)
EGM!(gl )
compute_distribution!(gl)
aggregate!(gl)

# 3) does the same as 2) in one step
compute_steady_state!(gl)

# 4) calibrate to target values for initial ss  -- requires 1)
calibrate!(gl)

# 5) calibrate to target values for terminal ss -- requires 3) or 4)
gl_tss = calibrate_terminal(gl)

 
  





# NOT WORKING YET: 

# re-calibrate (only phi2) to target for terminal ss
# EGM fails after ϕ is updated. Interpolation returns error meesage "knot-vectors must be sorted in increasing order"
gl_tss = gl
gl_tss.ϕ = gl_tss.ϕ₂
compute_tss!(gl_tss)

# compute transition dynamics
gl_trans = transition(gl,gl_tss)




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

 