# Replication of Guerrieri, Lorenzoni (QJE,2017)

> This replication study was part of our evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at SciencesPo Paris in Spring 2021

```@contents
```

In this replication study, we do ...


## Solving the Model

An easy way to solve the model in steady-state with a given parameterization is shown below. Note that all the steps have to occur sequentially.
First, instantiate the ModelGL structure that holds all parameter values and solutions of the model.
```
gl = ModelGL() 
```
Calling 
```
initialize!(gl)
``` 
On the object will set up the stationary distribution across productivity states and add the constraint explicitly to the bond grid.
Next,
```
EGM!(gl)
```
solves for the optimal policies.
```
compute_distribution!(gl)
```
calculates the joint distribution over productivity states and bond holdings and stores them in gl.JD.
```
aggregate!(gl)
```
computes aggregate values from the joint distribution.

 
## Calibrating the Model

To calibrate the model in steady-state to a given set of target values do the following:
First, instantiate the ModelGL structure that holds all parameter values, solutions and target values of the model.
```
gl = ModelGL() 
```
running the calibrate function will solve the steady-state model many times (steps from above) and update parameters in each iteration.
```
calibrate!(gl)
```
The solution calibrated parameters and solutions will be stored in the ModelGL structure.

## Functions

```@autodocs
Modules = [GL2017Replication]
```


end

