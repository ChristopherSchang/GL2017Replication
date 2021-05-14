# GL2017Replication

## About The Project

This replication study replicates parts of the paper `Credit crises, precautionary savings, and the liquidity trap` (Guerrieri, Lorenzoni (2017)). The code notation mostly follows the original version of the authors in MATLAB.

## Table of Contents
 
- [GL2017Replication](#gl2017replication)
  - [About The Project](#about-the-project)
  - [Table of Contents](#table-of-contents)
  - [Solving the Model](#solving-the-model)
  - [Calibrating the Steady-State](#calibrating-the-steady-state)
    - [Initial Steady-State](#initial-steady-state)
    - [Terminal Steady-State](#terminal-steady-state)
  - [Results](#results)
    - [Policies in Steady-State (initial)](#policies-in-steady-state-initial)
    - [Distribution in Steady-State (initial)](#distribution-in-steady-state-initial)
  - [Transition Dynamics](#transition-dynamics)
  - [References](#references)

## Solving the Model

An easy way to solve the model in steady-state with the given set of default parameters is shown below. Note that all the steps have to occur sequentially.
First, instantiate the ModelGL structure that holds all parameter values and solutions of the model.
```
gl = ModelGL() 
```
Calling 
```
compute_steady_state!(gl)
``` 
will solve for the steady-state policy function, joint distribution over productivity and asset states as well as aggregate variables.
 

## Calibrating the Steady-State 

### Initial Steady-State

To calibrate the model in steady-state to a given set of target values do the following:
First, instantiate the ModelGL structure as above.
```
gl = ModelGL() 
```
Calling the calibrate function will solve the model many times and calibrate the parameters to the target values:
```
calibrate!(gl)
```
The solution and the calibrated parameters will be stored in the ModelGL structure.

### Terminal Steady-State

To calibrate the model to the terminal steady-state we first need a calibrated model for the initial steady-state (step above):
```
gl = ModelGL()  
calibrate!(gl)
``` 
Given the solved initial steady-state we can calibrate the terminal steady-state to the new debt-target
```
gl_tss = calibrate_terminal(gl)
``` 
where `gl_tss` denotes the terminal steady-state object.

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

### Policies in Steady-State (initial)

<p float="left">
  <img src="images/c_pol_iss.png"    width="350"   />
  <img src="images/l_pol_iss.png"    width="350"   />  
</p> 

### Distribution in Steady-State (initial)

<p float="left"> 
  <img src="images/b_dist.png"      width="350" />
</p> 


## Transition Dynamics

## References 
Guerrieri, V. & Lorenzoni, G.
Credit crises, precautionary savings, and the liquidity trap 
The Quarterly Journal of Economics, Oxford University Press, 2017, 132, 1427-1467