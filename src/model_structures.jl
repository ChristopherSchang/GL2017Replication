
using Parameters
using QuantEcon
using DelimitedFiles
using MAT


# -----------------------------------------------------------------------------
# productivity process
# -----------------------------------------------------------------------------

nx      = 12
ρ       = 0.967  # persistence of productivity shock
σϵ      = sqrt(0.017)  # standard deviation  of productivity shock

#income process from MATLAB code
cd(@__DIR__)
pr = readdlm("ppr.csv",',',Float64)[:]
Pr = readdlm("Pr.csv",',',Float64) 
x  = readdlm("x.csv",',',Float64)[:]
steady_mat = matopen("steady.mat")
#read(steady_mat, "b_pol1")
 
θ               = [0; exp.(x)... ];
S               = length(θ);
fin             = 0.8820;        # job-finding probability
sep             = 0.0573;        # separation probability



"""
Structure holding all parameters and solutions to the steady-state model
"""
@with_kw mutable struct ModelGL

    # Deep parameters
    γ::Float64           = 4.       # risk aversion
    frisch::Float64      = 1.       # avg Frisch elast.
    r::Float64           = 2.5/400  # ss interest rate
    ρ::Float64           = ρ        # persistence of productivity shock
    σϵ::Float64          = σϵ       # standard deviation  of productivity shock  


    # Calibration targets
    NE::Float64         = 0.4       # avg hours of employed
    νY::Float64         = 0.4       # UI benefit per GDP
    B_4Y::Float64       = 1.6       # liquid wealth per annual GDP
    D1_4Y::Float64      = 0.18      # HH debt to annual GDP in initial ss
    D2_4Y::Float64      = 0.08      # HH debt to annual GDP in terminal ss

    η::Float64          = 1.5
    gameta::Float64     = γ / η;

    # Initial guesses for calibrated parameters, NE ~ Y
    β::Float64          = 0.8^0.25;                 # discount factor
    ν::Float64          = νY*NE;                    # UI benefits
    B::Float64          = B_4Y  * NE * 4;         # net supply of bonds
    ϕ₁::Float64         = D1_4Y * NE * 2;                  # borrowing constraint in initial ss
    ϕ₂::Float64         = D2_4Y * NE * 2;         # borrowing constraint in terminal ss
    ϕ::Float64          = ϕ₁
    ψ::Float64          = NE^(-γ) * (1-NE)^η;                  # disutility from labor as if representative agent

    # Numerical parameters
    maxit::Int64        = 500;   # maximum number of iterations in calibration
    cmin::Float64       = 1e-6;  # lower bound on consumption
    tol_pol ::Float64   = 1e-10; # tolerance lvl for policy function
    tol_dist::Float64   = 1e-10; # tolerance lvl for distribution
    tol_mkt ::Float64   = 1e-6;  # tolerance lvl for market clearing
    tol_cali::Float64   = 1e-4;  # tolerance lvl for calibration

    # bond (asset) grid
    nb::Int64                   = 200; # number of grid points
    Ic::Int64                   = 100;
    bmin::Float64               = -2.;   # lower bound
    bmax::Float64               = 50.;   # upper bound
    b_grid::Array{Float64,1}    = bmin .+ ((1:nb)/nb).^2 .* (bmax - bmin); # denser for low values
    db::Float64                 = 0.01; # step size for MPC

    # productivity shock process
    x::Array{Float64, 1}     = x
    #Pr_::Array{Float64, 2}   = Pr_
    #pr_::Array{Float64, 1}   = pr_

    # Add unemployment
    θ::Array{Float64, 1}= [0; exp.(x)... ];
    S::Int64            = length(θ);
    fin::Float64        = 0.8820;        # job-finding probability
    sep::Float64        = 0.0573;        # separation probability

    # new transition matrix
    Pr::Array{Float64, 2}   = Pr

    pr::Array{Float64,1}    = pr

    # Policies
    cl::Array{Float64,1}                = ones(S);     # consumption at borrowing constraint
    n_pol::Array{Float64,2}             = zeros(S, nb); # labor policy
    y_pol::Array{Float64,2}             = zeros(S, nb); # production policy
    b_pol::Array{Float64,2}             = zeros(S, nb); # savings policy
    c_pol::Array{Float64,2}             = zeros(S, nb); # consumption policy
    mpcs::Array{Float64,2}              = zeros(S, nb); # MPCs

    # Simulation
    JD::Array{Float64, 2}               = ones(S, nb) / (S+nb); # Joint compute_distribution

    # aggregates
    B_4Y_actual::Float64    = NaN
    D_4Y_actual::Float64    = NaN
    νY_actual::Float64      = NaN
    NE_actual::Float64      = NaN
    Y_actual::Float64       = NaN
    

    # diagnostic
    calibrated::Bool            = false
    steadystate_solved::Bool    = false
    steadystate_type::String    = "initial"
end


"""
Structure holding transition objects
"""
@with_kw mutable struct TransGL

    # SS Transition dynamics

    # Numerical parameters
    T::Int64               = 100                 # horizon, high enough for convergence to be complete
    maxit_trans::Int64     = 400                 # maximum number of iterations
    tol_mkt_trans::Float64 = 5e-6;               # tolerance for market clearing (lower than for ss computation)

    # Updating weights for interest rate
    speed::Float64 = 0.1
    decay::Float64 = 0.01
    weight_::Array{Float64,1} = exp.(-decay*(0:T-1))
    weight::Array{Float64,1}  = speed * weight_ / sum(weight_)

    # Transition objects
    S  = 13
    nb = 200
    ib_pol_t::Array{Int,3}     = zeros(S, nb, T); # sequence of policy functions
    wei_t::Array{Float64,3}    = zeros(S, nb, T); # sequence of weights on adjacent grid points
    Bdem_t::Array{Float64,1}   = zeros(T);        # bond demand
    Y_t::Array{Float64,1}      = zeros(T);        # GDP
    D_t::Array{Float64,1}      = zeros(T);        # debt
    D_4Y_t::Array{Float64,1}   = zeros(T);        # debt to GDP
    r_t::Array{Float64,1}      = zeros(T);        # interest rate path
    ϕ_t::Array{Float64,1}      = zeros(T+1);      # borrowing limit

    c_pol_t::Array{Float64,3} = zeros(S, nb, T);  # consumption policy in each period
    n_pol_t::Array{Float64,3} = zeros(S, nb, T);  # labour policy in each period
    y_pol_t::Array{Float64,3} = zeros(S, nb, T);  # income policy in each period

    JD_t::Array{Float64,3} = zeros(S,nb,T); # Joint compute_distribution at each point in time

end