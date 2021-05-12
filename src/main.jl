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

# save("gl.jld","gl",gl,"gl_tss",gl_tss)
# gl     = load("gl.jld")["gl"]
# gl_tss = load("gl.jld")["gl_tss"]

# 6) compute transition dynamics -- requires 4) and 5)

Tgl = TransGL()
gl_trans = transition!(gl,gl_tss,Tgl)


# -----------------------------------------------------------------------------
# plot results

using Plots
plotly(linewidth = 2.)
using Parameters, Statistics

@unpack b_grid,pr = gl                  # objects that are identical across both steady states
@unpack r_t,ϕ_t,Y_t,D_4Y_t = Tgl        # transition objects    

# Figure I (from paper)

Y1 = 0.4174 # !!! made up
plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.c_pol[2,b_grid .>= -gl.ϕ],label = "θ = 2")
plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.c_pol[8,b_grid .>= -gl.ϕ],
        linestyle = :dash,label = "θ = 8")
title!("consumption")
xaxis!([-gl.ϕ,12.5])
#savefig("images\\c_pol_iss.png")

plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.n_pol[2,gl.b_grid .>= -gl.ϕ],label =  "θ = 2")
plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.n_pol[8,gl.b_grid .>= -gl.ϕ],
        linestyle = :dash,label = "θ = 8")
title!("labor supply")
xaxis!([-gl.ϕ,12.5])
#savefig("images\\l_pol_iss.png")


# Figure IV (from paper)

b_acc1 =  ((pr'*gl.b_pol[:,b_grid     .>= -gl.ϕ])'     .- b_grid[b_grid .>= -gl.ϕ])/(4*Y1)
b_acc2 =  ((pr'*gl_tss.b_pol[:,b_grid .>= -gl_tss.ϕ])' .- b_grid[b_grid .>= -gl_tss.ϕ])/(4*Y1)
plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), b_acc1)
plot!(b_grid[b_grid .>= -gl_tss.ϕ]/(4*Y1), b_acc2,
        linestyle = :dash)
title!("bond accumulation policy")
xaxis!([-2,14],-2:2:14)
yaxis!([-0.4,0.6],-0.4:0.2:0.6)

bond_distribution1 = dropdims( pr'    *gl.JD[:,b_grid .>= -gl.ϕ]     ,dims=1)
bond_distribution2 = dropdims( pr'*gl_tss.JD[:,b_grid .>= -gl_tss.ϕ] ,dims=1)
plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1),  bond_distribution1)
plot!(b_grid[b_grid .>= -gl_tss.ϕ]/(4*Y1),  bond_distribution2,
        linestyle = :dash)
title!("bond distribution")
xaxis!([-2,14],-2:2:14)
yaxis!([0,0.004])
xaxis!([-ϕ,12.5])
#savefig("images\\b_dist_iss.png")

 # Figure III (from paper)
 Tp = 24    # number of periods plotted

 plot(0:Tp, ϕ_t[1:Tp+1]/(4*Y1))         # borrowing limit
 plot(0:Tp, D_4Y_t[1:Tp+1])             # debt2gdp ratio
 plot(0:Tp, [gl.r;r_t[1:Tp]]*400)       # annualized interest rate
 plot(0:Tp, [0, 100*(Y_t[1:Tp]/Y1-1)])   # output deviation from steady state


















 
