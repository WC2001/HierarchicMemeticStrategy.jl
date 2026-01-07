```@meta
CurrentModule = hms_jl
```

# hms_jl

The package **hms_jl** provides an implementation of the **Hierarchic Memetic Strategy (HMS)** — a population-based, multi-level metaheuristic for global optimization that combines evolutionary search with local optimization and hierarchical population management.

HMS organizes the search process into a **tree of demes**, where each level applies increasingly focused optimization strategies, enabling efficient exploration and exploitation of complex, multimodal objective functions.

## Getting started

To install the package, use:

```julia
] add hms_jl
```

Example of using the HMS algorithm to find minimum of the [Rosenbrock function](https://www.sfu.ca/~ssurjano/rosen.html) in 2 dimension version.


```@repl
using hms_jl

rosenbrock(x) = (1 - x[1])^2 + (100 * (x[2] - x[1]^2)^2)

problem = FunctionProblem(
    rosenbrock,
    [-30.0, -30.0],
    [30.0, 30.0],
    false
)

result = hms(optimization_problem = problem)
```
