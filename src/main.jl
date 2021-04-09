include("GL2017Replication.jl")
using .GL2017Replication



gl = ModelGL() 

# solve model
EGM!(gl)

# plot policy
using Plots
plotly()

plot(gl.b_grid,gl.b_pol[[1,5,12],:]',xlabel="bonds",ylabel="bond policy",label=["s=1" "s=5" "s=12"])
plot!(gl.b_grid,gl.b_grid,linestyle=:dot, label="45°")
