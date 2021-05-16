
using BenchmarkTools
using Parameters
using Roots
using Interpolations
using LinearAlgebra


function compute_steady_state!(gl::ModelGL) 
    initilize!(gl)
    EGM!(gl)
    compute_distribution!(gl)
    aggregate!(gl)
end


"""
Adds borrowing constraint as gridpoint
"""
function initilize!(gl::ModelGL)
 
    # Ensure that constraint is on the grid
    gl.b_grid = gl.bmin .+ ((1:gl.nb)/gl.nb).^2 .* (gl.bmax - gl.bmin); # denser for low values
    i = maximum(  [  searchsortedlast(gl.b_grid, -gl.ϕ) ,  1] )
    gl.b_grid[i] =  -gl.ϕ
end


"""
Calculates optimal policy functions with the EGM algorithm.
"""
function EGM!(gl::ModelGL)
    
    @unpack r,ϕ  = gl    
    @unpack β, ν,θ, η,ψ,γ, B,S,nb,gameta,fin,sep = gl               # model parameters
    @unpack cmin,tol_pol = gl                                       # numerical paras
    @unpack b_grid,Ic,db = gl                                       # grids
    @unpack pr,Pr = gl                                              # processes
    @unpack cl,n_pol,y_pol,b_pol,c_pol,tol_dist,mpcs = gl           # solutions

    # Update budget constraint
    τ   = (pr[1]*ν + r/(1+r)*B) / (1 - pr[1]);      # labor tax
    z   = [ν, -τ.*ones(S-1)...];                    # full transfer scheme (tau tilde in paper)

    # Find consumption at the lower bound of the state space
    fac = (ψ ./ θ) .^ (1/η);
    for s = 1:S
        
        cl[s] =   find_zero( x -> find_cl(x,s, -ϕ, -ϕ, r, θ, z, fac, gameta)
                                , (cmin, 100), Bisection()  )
    end

    # A) Solve for consumption policy
    #----------------------------------
    # Initial guess for policy function
    c_pol[:,:]  = r .* ones(S,nb) .* b_grid'
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
        itp = extrapolate(    interpolate((b,),c, Gridded(Linear())), Interpolations.Flat() )
        c_poli[s,:] = itp(  b_grid   )
        # interp1(b, c, b_grid, 'linear', 'extrap');
      end

      # check convergence
      c_poli[c_poli .<= cmin] .= cmin
      dif    = norm(c_poli - c_pol)/(1+norm(c_pol) );

      # update
      c_pol[:,:] = copy(c_poli);
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

""" Computes the joint distribution of productivity and bond holdings.
Requires that optimal policies have already been solved (EGM!)
"""
function compute_distribution!(gl::ModelGL)

    @unpack b_pol, b_grid,tol_dist,S,nb, JD,Pr,fin,pr,sep,ϕ = gl


    # B) Find invariant distribution
    #----------------------------------
    # Assign weights to adjacent grid points proportionally to distance
    idxes       = [ searchsortedlast(b_grid,b_pol[s, b]) for s = 1:size(b_pol)[1], b = 1:size(b_pol)[2] ]
    weights     = zeros(size(b_pol))
    for s = 1:size(b_pol)[1]
    for b = 1:size(b_pol)[2]
        idx             = idxes[s,b]
        idx_prime       = minimum(  [ idx+1,size(b_pol)[2] ]    )
        #db              = maximum(  [ b_grid[idx_prime] - b_grid[idx]  , eps() ] )
        db              =  b_grid[idx_prime] - b_grid[idx]
        weights[s,b]    = ( b_pol[s,b]  - b_grid[idx] ) / db;
    end
    end

    # Iterate asset transition matrix starting from uniform distribution
    dif = 1;
    iter =1
    JD = ones(S,nb) ./ (S*nb);
    while (dif > tol_dist) && iter < 1000

        iter += 1
        JDp = zeros( size(JD) )
        for s = 1:S
        for b = 1:nb
            for si = 1:S
                JDp[si, idxes[s, b]]     = (1 - weights[s, b]) * Pr[s, si] * JD[s, b]  + JDp[si, idxes[s, b]]
                JDp[si, idxes[s, b] + 1] =      weights[s, b]  * Pr[s, si] * JD[s, b]  + JDp[si, idxes[s, b] + 1]
            end
        end
        end

        # check convergence
        dif = norm( JDp - JD ) /( 1+norm(JD) )

        # make sure that distribution integrates to 1
        JD[:,:]  = JDp / sum( JDp );

    end

    gl.JD = copy(JD)
    gl.steadystate_solved = true
    return JD
end

"""
Computes aggregate values and their distance from target.
Requires that the joint distribution has already been solved (computer_distribution!)
"""
function aggregate!(gl::ModelGL;terminal::Bool = false)

    # Compute Steady-state aggregates
    @unpack ν,β,ϕ,ψ,η  = gl
    @unpack NE,νY,B_4Y,D1_4Y,D2_4Y  = gl
    @unpack JD,b_grid,B,y_pol,n_pol  = gl
    @unpack tol_cali,tol_mkt  = gl

    # C) Check market clearing and calibration
    #-------------------------------------------------
    # Bond market clearing, i refers to current iteration
    Bi      = (sum(JD,dims=1) * b_grid)[1]
    res_mkt = abs(B - Bi)

    # Calibration statistics, i refers to current iteration
    Yi    = sum(JD  .* y_pol  )                                     # GDP
    NEi   = sum(  JD .* (n_pol.*(n_pol.>0) ) )  / sum(  JD .*(n_pol.>0)  )       # avg hours of employed
    B_4Yi =  Bi / Yi / 4;                                                # debt ratio
    Di    = - (   sum(JD,dims=1) *  ( b_grid .* (b_grid .< 0.) )  )[1]
    D_4Yi =  Di / Yi / 4;                                                # debt ratio
    νYi    =  ν / Yi;                                                    # UI benefit ratio
    if terminal == false
        res_cali = maximum(abs.([B_4Yi, D_4Yi, νYi, NEi] .- [B_4Y, D1_4Y, νY, NE]));
    else
        res_cali = maximum(abs.([  D_4Yi ] .- [  D2_4Y ]));
    end

    #store actual values 
    gl.B_4Y_actual  = B_4Yi
    gl.D_4Y_actual  = D_4Yi
    gl.νY_actual    = νYi
    gl.NE_actual    = NEi
    gl.Y_actual     = Yi
    return Bi ,res_mkt ,Yi ,NEi,B_4Yi ,Di,D_4Yi,νYi,res_cali
end

"""
calibrate the model to the target values.
Does not require any prior function calls.
"""
function calibrate!(gl::ModelGL)

    @unpack NE,νY,η,B_4Y,D1_4Y  = gl
    @unpack tol_cali,tol_mkt  = gl
    ite = 0
    while true
        ite +=1
        initilize!(gl)
        EGM!(gl)
        compute_distribution!(gl)
        Bi ,res_mkt ,Yi ,NEi,B_4Yi ,Di,D_4Yi,νYi,res_cali = aggregate!(gl)
        println( string(ite)*"  "*string(res_mkt)*"   "*string(res_cali) )
        # Check convergence of both, update if necessary
        if (res_cali < tol_cali) && (res_mkt < tol_mkt)
            break
        else
            # discount factor
            btemp   = -log(1- gl.β);
            btemp   = btemp - .1*(Bi - B_4Y*4*Yi);
            gl.β       = 1 - exp(-btemp);

            # rest will be updated based on these
            ϕ_d    = gl.ϕ * D1_4Y / D_4Yi;
            ν_d    = νY * Yi;
            ψ_d    = gl.ψ * ((1-NE) / (1-NEi))^η;

            gl.ϕ   = gl.ϕ  + 0.1 * (ϕ_d   - gl.ϕ);
            gl.ν   = gl.ν  + 1   * (ν_d   - gl.ν);
            gl.ψ   = gl.ψ  + 1   * (ψ_d   - gl.ψ);

            # Update aggregates
            gl.B = gl.B + 0.1 * (Bi - gl.B);
        end
    end
    gl.calibrated = true
end


"""
calibrate the model to the terminal target value (only borrowing constraint)
Requires a calibrated initital steady-state object.
Returns a new terminal steady-state object.
"""
function calibrate_terminal(gl_initial::ModelGL)
    
    gl = deepcopy(gl_initial)   # terminal steadystate
    gl.ϕ = gl.ϕ₂                # update initial guess

    @unpack D2_4Y  = gl
    @unpack tol_cali,tol_mkt  = gl
    ite = 0
    while true
        ite +=1
        initilize!(gl)
        EGM!(gl)
        compute_distribution!(gl)
        Bi ,res_mkt ,Yi ,NEi,B_4Yi ,Di,D_4Yi,νYi,res_cali = aggregate!(gl,terminal = true)
        println( string(ite)*"  "*string(res_mkt)*"   "*string(res_cali) )
        # Check convergence of both, update if necessary
        if (res_cali < tol_cali) && (res_mkt < tol_mkt)
            break
        else 
            gl.r = gl.r - 1/400*(Bi-gl.B)
            
            ϕ_d    = gl.ϕ * D2_4Y / D_4Yi;  
            gl.ϕ   = gl.ϕ  + 0.1 * (ϕ_d   - gl.ϕ); 
 
        end
    end
    gl.calibrated = true
    gl.steadystate_type = "terminal"
    return gl
end


"""
compute transition path.
"""

function transition!(gl::ModelGL,gl_tss::ModelGL,Tgl::TransGL)
    @unpack T, maxit_trans, tol_mkt_trans, weight, JD_t = Tgl
    @unpack Pr,B,nb = gl        

    gl_t = deepcopy(gl)  # define new object for each time period

    # Define credit crunch parameters
    t_shock = 6                           # quarters to get to new constraint
    dϕ      = (gl.ϕ - gl_tss.ϕ)/t_shock   # step size
    Tgl.ϕ_t  = [max(gl_tss.ϕ,gl.ϕ .- dϕ*t) for t in 0:T] # full path

    # Initial guess of r_t
    Tgl.r_t = gl_tss.r .* ones(T)  # interest rate in terminal ss
    for it = 1:maxit_trans
    
        # 1) Iterate HH problem backwards      
        gl_t.c_pol = gl_tss.c_pol   # start from terminal ss
        
        for t = T:-1:1
            # current value of r and ϕ
            gl_t.r = Tgl.r_t[t]
            gl_t.ϕ = Tgl.ϕ_t[t+1]
            # Ensure that constraint is on the grid
            initilize!(gl_t)
            # Compute policy functions and save distribution weights
            EGM_trans!(gl_t,Tgl,t)
        end

        # 2) Iterate distribution forward

        JD = gl.JD;  # start from initial ss
        for t = 1:T
            # store distribution over assets and productivity states
            JD_t[:,:,t] .= JD

            # Ensure that constraint is on the grid
            gl_t.ϕ = Tgl.ϕ_t[t+1]
            initilize!(gl_t)

            # calculate bond market clearing
            Tgl.Bdem_t[t] = (sum(JD,dims=1) * gl_t.b_grid)[1]

            # calculate other aggregates
            Tgl.Y_t[t]    = sum(JD  .* Tgl.y_pol_t[:,:,t])                         # GDP
            Tgl.D_t[t]    = -(sum(JD,dims=1) * (gl_t.b_grid .* (gl_t.b_grid .< 0.)))[1]  # Debt
            Tgl.D_4Y_t[t] =  Tgl.D_t[t] / Tgl.Y_t[t] / 4                               # Debt2GDP annual           

            # iterate on distribution
            JDp = zeros(size(JD))
            for s = 1:S
            for b = 1:nb
                for si = 1:S
                    JDp[si, Tgl.ib_pol_t[s, b, t]]     = (1 - Tgl.wei_t[s, b, t]) * Pr[s, si] * JD[s, b]  + JDp[si, Tgl.ib_pol_t[s, b, t]]
                    JDp[si, Tgl.ib_pol_t[s, b, t] + 1] =      Tgl.wei_t[s, b, t]  * Pr[s, si] * JD[s, b]  + JDp[si, Tgl.ib_pol_t[s, b, t] + 1]
                end
            end
            end

            # make sure that distribution integrates to 1
            JD[:,:]  = JDp / sum(JDp)            

        end
        
        # 3) Check convergence
        
        # Market clearing
        BM_t = Tgl.Bdem_t .- B

        # Metric of deviation
        res_BM = sqrt(BM_t' * BM_t / T)

        # Report progress
        println( string(it)*"  "*string(res_BM))

        if res_BM < tol_mkt_trans
            break
        end

        # Update 
        Tgl.r_t = Tgl.r_t - weight.*BM_t
    end

end


"""
Calculates optimal policy functions with one step of the EGM algorithm.
"""
function EGM_trans!(gl::ModelGL,Tgl::TransGL,t::Int)
    
        @unpack r,ϕ  = gl    
        @unpack β, ν,θ, η,ψ,γ, B,S,gameta = gl            # model parameters
        @unpack cmin = gl                                 # numerical paras
        @unpack b_grid,Ic,bmax = gl                       # grids
        @unpack pr,Pr = gl                                # processes
        @unpack cl, c_pol = gl                            # solutions

        # Update budget constraint
        τ   = (pr[1]*ν + r/(1+r)*B) / (1 - pr[1]);      # labor tax
        z   = [ν, -τ.*ones(S-1)...];                    # full transfer scheme (tau tilde in paper)

        # Find consumption at the lower bound of the state space
        fac = (ψ ./ θ) .^ (1/η);
        for s = 1:S
                cl[s] =   find_zero( x -> find_cl(x,s, -Tgl.ϕ_t[t], -Tgl.ϕ_t[t+1], r, θ, z, fac, gameta)
                , (cmin, 100), Bisection()  )
        end

        # expected marginal utility tomorrow
        ui = Pr *  c_pol.^(-γ) ;
        ui = ui[:, b_grid .>= -ϕ];

        # previous consumption and savings policy
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

            # update policies
            itp = extrapolate(    interpolate((b,),c, Gridded(Linear())), Interpolations.Flat() )
            gl.c_pol[s,:] = itp(b_grid)
            gl.c_pol[c_pol .<= cmin] .= cmin

            npp   = 1 .- fac[s] .*gl.c_pol[s,:] .^gameta
            npp[npp .<= 0.] .= 0.
            gl.n_pol[s,:] .= npp
            gl.y_pol[s,:] = θ[s] .* gl.n_pol[s,:];
            bpp        = (1+r) * (b_grid .+ gl.y_pol[s,:] .- gl.c_pol[s,:] .+ z[s])
            bpp[bpp .<= -ϕ] .= -ϕ
            bpp[bpp .>= bmax] .= bmax - 0.000001
            gl.b_pol[s,:] .= bpp

            # save policies
            Tgl.c_pol_t[s, :, t] = gl.c_pol[s, :];
            Tgl.n_pol_t[s, :, t] = gl.n_pol[s, :];
            Tgl.y_pol_t[s, :, t] = gl.y_pol[s, :];
        end
        
        # save weights needed to iterate distribution
        idxes       = [searchsortedlast(b_grid,gl.b_pol[s, b]) for s = 1:size(gl.b_pol)[1], b = 1:size(gl.b_pol)[2] ]
        weights     = zeros(size(gl.b_pol))
        for s = 1:size(gl.b_pol)[1]
        for b = 1:size(gl.b_pol)[2]
            idx             = idxes[s,b]
            idx_prime       = minimum(  [ idx+1,size(gl.b_pol)[2] ]    )
            db              =  b_grid[idx_prime] - b_grid[idx]
            weights[s,b]    = ( gl.b_pol[s,b]  - b_grid[idx] ) / db;
        end
        end
        Tgl.ib_pol_t[:,:,t] = idxes
        Tgl.wei_t[:,:,t]    = weights 
end

"""
Helper function to find consumption at the constraint.
Returns the residual the budget constraint
"""
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





