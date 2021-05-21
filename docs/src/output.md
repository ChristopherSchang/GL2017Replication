# Output

## Results
At any point we can query basic statistics (parameter values, aggregates and targets) with the describe function
```
describe(gl)
``` 
or
```
describe(gl,gl_tss)
``` 
for a side-by-side comparison of the initial and terminal steady-state.

## Plotting
To plot some of the key figures you can use the following plotting functions
```
plots_figure_1(gl)
plots_figure_3(gl,Tgl)
plots_figure_4(gl,gl_tss)
```
for Figures I, III, and IV, respectively. The arguments are again the initial and terminal steady-state object (gl,gl_tss) and the transition object (Tgl). 

