```@meta
CurrentModule = HierarchicMemeticStrategy
```

# HierarchicMemeticStrategy.jl

The package **HierarchicMemeticStrategy.jl** provides an implementation of the **Hierarchic Memetic Strategy (HMS)** — a population-based, multi-level metaheuristic for global optimization that combines evolutionary search with local optimization and hierarchical population management.

HMS organizes the search process into a **tree of demes**, where each level applies increasingly focused optimization strategies, enabling efficient exploration and exploitation of complex, multimodal objective functions.

## Getting started

To install the package, use:

```julia
] add HierarchicMemeticStrategy
```

Example of using the HMS algorithm to find minimum of the [Rosenbrock function](https://www.sfu.ca/~ssurjano/rosen.html) in 2 dimension version.


```@repl
using HierarchicMemeticStrategy

rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)

problem = FunctionProblem(
    fitness_function=rosenbrock,
    lower=[-30.0, -30.0],
    upper=[30.0, 30.0]
)

result = hms(optimization_problem=problem)
```
