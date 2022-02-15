# Thyrosim

`Thyrosim.jl` is an updated version of [THYROSIM](http://biocyb1.cs.ucla.edu/thyrosim/cgi-bin/Thyrosim.cgi). It heavily relies on the amazing [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package to solve systems of ODEs. 

This tool is developed to

1. Getter optimize replacement LT4 and LT4+LT3 dosing for treating hypothyroid patients, based on their individual hormone levels, heights, weights and gender
2. Aid in understanding more precisely how gender and BMI impact thyroid hormone transient and steady state dynamical regulation over time in these patients.

## Installation

Download and install [Julia](https://julialang.org/downloads/). This package supports Julia `v1.6`+ for Mac, Linux, and window machines. Within Julia, copy and paste the following:
```julia
using Pkg
pkg"add https://github.com/biona001/Thyrosim.jl"
```
For running the examples below, the following packages are also necessary
```julia
pkg"add Plots"
```

## Example usages

Load packages
```julia
using Thyrosim, Plots
```

Healthy female patient simulation
```julia
sex = false # true = male, false = female
h = 1.63    # unit = meters
w = 60      # unit = KG
sol = simulate(h, w, sex, days=50, warmup=true) # solve ODE
plt = output_plot(sol, title="Female, $h meters, $w KG") # make plot
```
![sim1](images/sim1.png)

## Citation and Reproducibility

Coming soon

## Bug reports 

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error.

PRs and feature requests are welcomed!
