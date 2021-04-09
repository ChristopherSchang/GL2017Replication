
using Parameters
using QuantEcon
 

 
# -----------------------------------------------------------------------------
# productivity process
# -----------------------------------------------------------------------------

nx      = 12
ρ       = 0.967  # persistence of productivity shock
σϵ      = 0.017  # standard deviation  of productivity shock

mc_productivity  = tauchen(nx, ρ, σϵ, 0.0, 3)
x               = mc_productivity.state_values
Pr              = mc_productivity.p
pr            = stationary_distributions(mc_productivity)[1]

# !!! cannot replicate exact values !!!
# -----------------------------------------------------------------------------


""" 
single structure holding all parameters and solutions to the model 
"""
@with_kw mutable struct ModelGL

    # Deep parameters
    γ::Float64           = 4.       # risk aversion
    frisch::Float64      = 1.       # avg Frisch elast.
    r::Float64           = 2.5/400  # ss interest rate  
    

    # Calibration targets
    NE::Float64         = 0.4       # avg hours of employed
    νY::Float64         = 0.4       # UI benefit per GDP
    B_4Y::Float64       = 1.6       # liquid wealth per annual GDP
    D1_4Y::Float64      = 0.18      # HH debt to annual GDP in initial ss
    D2_4Y::Float64      = 0.08      # HH debt to annual GDP in terminal ss

    η::Float64          = 1/frisch * (1 - NE) / NE
    gameta::Float64     = γ / η;

    # Initial guesses for calibrated parameters, NE ~ Y
    β::Float64          = 0.8^(1/4);              # discount factor
    ν::Float64          = νY    * NE;             # UI benefits
    B::Float64          = B_4Y  * NE * 4;         # net supply of bonds
    ϕ₁::Float64         = D1_4Y * NE * 2;         # borrowing constraint in initial ss
    ϕ₂::Float64         = D2_4Y * NE * 2;         # borrowing constraint in terminal ss
    pssi::Float64       = NE^(-γ) * (1-NE)^η;     # disutility from labor as if representative agent

    # Numerical parameters
    maxit::Int64        = 500;   # maximum number of iterations in calibration
    cmin::Float64       = 1e-6;  # lower bound on consumption
    tol_pol ::Float64   = 1e-10; # tolerance lvl for policy function
    tol_dist::Float64   = 1e-10; # tolerance lvl for distribution
    tol_mkt ::Float64   = 1e-6;  # tolerance lvl for market clearing
    tol_cali::Float64   = 1e-4;  # tolerance lvl for calibration
 
    
    # bond (asset) grid
    nb::Int64                   = 200; # number of grid points
    Ic                          = 100;
    bmin                        = -2;   # lower bound
    bmax                        = 50;   # upper bound
    b_grid::Array{Float64,1}    = bmin .+ ((1:nb)/nb).^2 .* (bmax - bmin); # denser for low values
    db                          = 0.01; # step size for MPC

    # productivity shock process
    x::Array{Float64, 1}    = x
    Pr::Array{Float64, 2}   = Pr
    pr::Array{Float64, 1}   = pr
    
    # Add unemployment
    θ       = [0; exp.(x)... ];
    S       = length(θ);
    fin     = 0.8820;        # job-finding probability
    sep     = 0.0573;        # separation probability


    # Policies
    cl::Array{Float64,1}                = ones(S);     # consumption at borrowing constraint
    n_pol::Array{Float64,2}             = zeros(S, nb); # labor policy
    y_pol::Array{Float64,2}             = zeros(S, nb); # production policy
    b_pol::Array{Float64,2}             = zeros(S, nb); # savings policy
    c_pol::Array{Float64,2}             = zeros(S, nb); # consumption policy
    mpcs::Array{Float64,2}              = zeros(S, nb); # MPCs

    # Simulation
    JD::Array{Float64, 2}               = ones(S, nb) / (S+nb); # Joint compute_distribution
    

 
end

 
 