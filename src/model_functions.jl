
using BenchmarkTools 
using Parameters
using Roots
using Interpolations
using LinearAlgebra

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
    
    @unpack β, ν,θ, ϕ₁, ϕ₂,η,pssi,γ, B,S,r,nb,gameta,fin,sep = gl    # model parameters
    @unpack cmin,tol_pol = gl                           # numerical paras
    @unpack b_grid,Ic,db = gl                                 # grids
    
    @unpack pr,Pr = gl                                     # processes
    @unpack cl,n_pol,y_pol,b_pol,tol_dist,mpcs = gl                   # solutions
    
    
    # new transition matrix
    Pr = [1-fin  fin.*pr'; 
            sep .*ones(S-1)  (1-sep) .*Pr];
    
    # find new invariate distribution
    pr  = [0, pr...];
    dif = 1;
    while dif > tol_dist
        pri = (pr'*Pr)';
        dif = maximum( abs.(pri .- pr) );
        pr  = pri;
    end
    
    ϕ =  ϕ₁
        
    # Update budget constraint
    τ   = (pr[1]*ν + r/(1+r)*B) / (1 - pr[1]);  # labor tax
    z   = [ν, -τ.*ones(S-1)...];                   # full transfer scheme (tau tilde in paper)
    
    # Find consumption at the lower bound of the state space
    fac = (pssi ./ θ) .^ (1/η); 
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
    
    # policy Convergence
    dif = 1;
    while dif > tol_pol
      # expected marginal utility tomorrow
      ui = Pr *  c_pol.^(-γ) ;
      ui = ui[:, b_grid .>= -ϕ];
    
      for s=1:S
        # unconstrained
        c = ((1+r) * β * ui[s,:]) .^ (-1/γ);                   # Euler
        n =   (1 .- fac[s]*c.^gameta )                         # labor supply
        n = [ n[i] < 0. ? 0. : n[i] for i = 1:length(n)]
        b = b_grid[b_grid .>= -ϕ] / (1+r) .+ c .- θ[s] .*n .- z[s]; # budget
    
        # constrained
        if b[1] > -ϕ
            c_c = range(cl[s], stop=c[1], length=Ic);
            n_c =   (1 .- fac[s]*c_c.^gameta )                         # labor supply
            n_c = [ n_c[i] < 0. ? 0. : n_c[i] for i = 1:length(n_c)]    
            b_c = -ϕ/(1+r) .+ c_c .- θ[s] .*n_c .- z[s]; # budget
            b   = [b_c[1:Ic-1]..., b...];
            c   = [c_c[1:Ic-1]..., c...];
        end
        # !!! fix order !!!
        sort!(b)
        itp = extrapolate(    interpolate((b,),c, Gridded(Linear())), Interpolations.Flat() )
        c_poli[s,:] = itp(  b_grid   )
        # interp1(b, c, b_grid, 'linear', 'extrap');
      end
    
      # check convergence
      c_poli[c_poli .<= cmin] .= cmin
      dif    = norm(c_poli - c_pol)/(1+norm(c_pol) );
    
      # update
      c_pol = copy(c_poli);
    end
    
    # Save other policy functions and MPCs
    for s = 1:S
      npp   = 1 .- fac[s] .*c_pol[s,:] .^gameta
      npp[npp .<= 0.] .= 0.
      n_pol[s,:] .= npp
      y_pol[s,:] = θ[s] .* n_pol[s,:];
      bpp        = (1+r) * (b_grid .+ y_pol[s,:] .- c_pol[s,:] .+ z[s])
      bpp[bpp .<= -ϕ] .= -ϕ
      b_pol[s,:] .= bpp
      itp = extrapolate(    interpolate((b_grid,),c_pol[s,:], Gridded(Linear())), Interpolations.Flat() )
      mpcs[s,:]  = (    itp( b_grid .+ db )  .- c_pol[s,:]  ) ./ db;
      # (interp1(b_grid, c_pol(s,:), b_grid + db, 'linear', 'extrap') - c_pol(s,:)) / db;
    end


end



function find_cl(c, j, b1, b2, r, θ, z, fac, gameta)
    # Find consumption for current assets b1 and next assets b2 

    n = max(0, 1 - fac[j]*c.^gameta);
    F = b1 - b2/(1+r) - c + θ[j]*n + z[j];
    
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
