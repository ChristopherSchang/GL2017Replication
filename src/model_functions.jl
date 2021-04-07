

using BenchmarkTools
using Roots



""" Function to solve steady-state objects.  
Inputs
    p::ModelParameters          :   Structure holding parameter Values
    grid::ModelSolutions        :   Structure for grids
    solution::ModelSolutions    :   Structure for solutions
"""
function solve_steady_state!(gl::ModelGL)
    #do stuff

end

function EGM!(gl::ModelGL)
    
    @unpack ν,θ,ϕ,η,pssi, B,S,r,nb = gl
    @unpack cmin = gl
    @unpack pr = gl
    @unpack cl = gl

    # Update budget constraint
    τ   = (pr[1]*ν + r/(1+r)*B) / (1 - pr[1]);  # labor tax
    z   = [ν, -τ.*ones(S-1)];                   # full transfer scheme (tau tilde in paper)

    # Find consumption at the lower bound of the state space
    fac = (pssi / θ) ^ (1/η); 
    for s = 1:S 
        cl[s] =   find_zero( x -> find_cl(x,s, -ϕ, -ϕ, r, θ, z, fac, gameta) 
                                , (cmin, 100), Bisection()  )
    end

    # A) Solve for consumption policy
    #----------------------------------
    # Initial guess for policy function
    c_pol  = r .* ones(S,nb) .* b_grid'  
    c_pol[c_pol .<= cmin] .= cmin
    c_poli = copy(c_pol)                         

    # .... to be continued....

end



function find_cl(c, j, b1, b2, r, θ, z, fac, gameta)
    # Find consumption for current assets b1 and next assets b2 

    n = max(0, 1 - fac(j)*c.^gameta);
    F = b1 - b2/(1+r) - c + θ[j]*n + z(j);
    
    return F
end




"""
    hello(who::String)

Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)

Return `x + 5`.
"""
domath(x::Number) = x + 5
