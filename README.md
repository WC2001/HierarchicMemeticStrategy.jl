# HierarchicMemeticStrategy.jl

<!-- badges: start -->

<!-- badges: end -->

The HMS (Hierarchic Memetic Strategy) is a global optimization strategy consisting of a multi-population evolutionary strategy. The HMS stores a hierarchic data structure of component populations
as a tree with a fixed maximal height and a variable internal node degree. Each component population is processed by a dedicated evolutionary algorithm. This package offers a concise Julia implementation of HMS, including examples that showcase different genetic algorithms usage.

### Literature

-   J. Sawicki, M. Łoś, M. Smołka, R. Schaefer. Understanding measure-driven algorithms solving irreversibly ill-conditioned problems. Natural Computing 21:289-315, 2022. doi: [10.1007/s11047-020-09836-w](https://doi.org/10.1007/s11047-020-09836-w)

## Installation

To install the package, use:

```julia
] add HierarchicMemeticStrategy
```

## Usage

To run the HMS with a default configuration provide `FunctionProblem` 
with fitness function and bounds.

``` julia
using HierarchicMemeticStrategy

rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)

problem = FunctionProblem(
    fitness_function=rosenbrock,
    lower=[-30.0, -30.0],
    upper=[30.0, 30.0]
)

result = hms(optimization_problem=problem)
```

### Adjusted configuration

```julia
# optimizing eggholder function
function eggholder(x::Vector{Float64})
    x1, x2 = x[1], x[2]
    return -(x2 + 47) * sin(sqrt(abs(x2 + x1/2 + 47))) - x1 * sin(sqrt(abs(x1 - (x2 + 47))))
end

lower = [-512.0, -512.0]
upper = [512.0, 512.0]
problem = FunctionProblem(
    fitness_function=eggholder,
    lower=lower,
    upper=upper
)

seed = 42
sigma = [[100.0, 100.0], [60.0, 60.0]]
level_config = [
    TreeLevelConfig(EvolutionaryGAMetaepoch, Dict("seed" => seed)),
    TreeLevelConfig(EvolutionaryCMAESMetaepoch, Dict("seed" => seed)),
]

result = hms(
    optimization_problem=problem,
    level_config=level_config,
    sigma=sigma,
    seed=seed 
)
```
