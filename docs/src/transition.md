# Transition

This part is of particular interest. When we tried to replicate Figure 3 of the original paper without loading the transition objects previously obtained by the authors (i.e., we set $rerun \neq 1$), but keeping the default values of $speed = 0.5$ and $decay = 0.3$ (two parameters governing the updating rule of ), we couldn't make the code converge. It seems to us that the authors used two different parameters that are a priori impossible to deduce. Importantly, for one of the few combinations of $speed$ and $decay$ that would produce output (i.e., 0.1 and 0.01) in the original MATLAB file, we obtain essentially the same evolution for all variables included in Figure 3 (trivially, the borrowing constraint, but also household debt, the interest rate, and aggregate output). This somehow reassures us about the validity of our approach, although we cannot obviously rule out completely the possibility that we missed some other features of the code explaining the puzzle.

end

