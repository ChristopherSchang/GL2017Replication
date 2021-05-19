# Replication of Guerrieri, Lorenzoni (QJE,2017)

> This replication study was part of our evaluation for the course [Numerical Methods](https://floswald.github.io/NumericalMethods/) at SciencesPo Paris in Spring 2021

```@contents
```

## About The Project

This replication study replicates parts of the paper ´Credit crises, precautionary savings, and the liquidity trap´ (Guerrieri, Lorenzoni (2017)). The code notation mostly follows the original version of the authors in MATLAB.

## Structure
The Package is build around two structures: 
1. The ModelGL structure that holds all parameters and (if solved/calibrated) the policy functions and aggregate variables.
```
gl = ModelGL() 
```
At any point you query the status of the model with
```
describe(gl)
```
2. The TransGL structure that holds the objects for the transition from an initital to a terminal steady-state
```
gl_trans = TransGL() 
```

## References 
Guerrieri, V. & Lorenzoni, G.
Credit crises, precautionary savings, and the liquidity trap 
The Quarterly Journal of Economics, Oxford University Press, 2017, 132, 1427-1467

end

