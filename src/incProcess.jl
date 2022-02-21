"""
# new transition matrix
mc_productivity  = tauchen(nx, ρ, σϵ, 0.0, 2)
x               = mc_productivity.state_values
Pr_             = mc_productivity.p
pr_              = stationary_distributions(mc_productivity)[1]

Pr              = [1-fin              fin.*pr_';
                    sep .*ones(S-1)  (1-sep) .*Pr_];
    
function initt(Pr,pr_)
    tol_dist    = 1e-10;
    # find new invariate distribution
    pr  = [0, pr_...];
    dif = 1;
    while dif > tol_dist 
        pri = (pr'*Pr)';
        dif = maximum( abs.(pri .-  pr) );
        pr  = pri;
    end
    return pr
end

pr = initt(Pr,pr_)
 
"""