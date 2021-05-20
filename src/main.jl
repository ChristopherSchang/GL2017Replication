using Revise, JLD
include("GL2017Replication.jl")
using .GL2017Replication

# -----------------------------------------------------------------------------

# 1) Load model structure
gl = ModelGL()
describe(gl)
 
# 2) solve model with default parameters  -- requires 1)
compute_steady_state!(gl)
describe(gl)

# 3) calibrate to target values for initial ss  -- requires 1)
calibrate!(gl)
describe(gl)

# 4) calibrate to target values for terminal ss -- requires 3)  
gl_tss = calibrate_terminal(gl)
describe(gl,gl_tss)

# save("gl.jld","gl",gl,"gl_tss",gl_tss)
# gl     = load("gl.jld")["gl"]
# gl_tss = load("gl.jld")["gl_tss"]

# 5) compute transition dynamics -- requires 3) and 4)

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

Y1 = gl.Y_actual
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
#savefig("images\\b_acc.png")

bond_distribution1 = dropdims( pr'    *gl.JD[:,b_grid .>= -gl.ϕ]     ,dims=1)
bond_distribution2 = dropdims( pr'*gl_tss.JD[:,b_grid .>= -gl_tss.ϕ] ,dims=1)
plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1),  bond_distribution1,label = "initial")
plot!(b_grid[b_grid .>= -gl_tss.ϕ]/(4*Y1),  bond_distribution2,label = "terminal",
        linestyle = :dash)
title!("bond distribution")
xaxis!([-2,14],-2:2:14)
yaxis!([0,0.004])
xaxis!([-gl.ϕ,12.5])
#savefig("images\\b_dist.png")

 # Figure III (from paper)
 Tp = 24    # number of periods plotted

 p1 = plot(0:Tp, Tgl.ϕ_t[1:Tp+1]./(4*Tgl.Y_t[1:Tp+1]),ylims = (0.5,1), yticks = 0.5:0.1:10, legend = false, title = "borrowing limit")         # borrowing limit
 p2 = plot(0:Tp, Tgl.D_4Y_t[1:Tp+1], ylims = (0.08,0.2), yticks = 0:0.02:0.2, legend = false, title = "household debt-to-GDP ratio")             # debt2gdp ratio
 p3 = plot(0:Tp, [gl.r;Tgl.r_t[1:Tp]].*400, ylims = (-2,2.5), yticks = -2:0.5:2, legend = false, title = "interest rate")       # annualized interest rate
 p4 = plot(0:Tp, [0, 100*(Tgl.Y_t[1:Tp+1]./Y1.-1)], ylims = (-1.2,0), yticks = -1:0.2:0, legend = false, title = "output")   # output deviation from steady state
 plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
#savefig("images\\trans.png")
















 
