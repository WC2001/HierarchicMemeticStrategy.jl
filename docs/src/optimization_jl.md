# Optimization.jl Interface

`HierarchicMemeticStrategy.jl` is integrated with the [Optimization.jl](https://docs.sciml.ai/Optimization/stable/) (SciML) ecosystem. This allows you to utilize a unified syntax for defining optimization problems and provides an easy way to swap between HMS and other global or local solvers.

## Basic Usage

To use HMS within this ecosystem, you define the `OptimizationProblem` and pass `HMSSolver()` as the solver to the `solve` function.

### Example: Solving the Rosenbrock Function

```@repl optim
using Optimization
using HierarchicMemeticStrategy

# 1. Define the objective function (Optimization.jl expects f(x, p))
rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2;

# 2. Define starting point, parameters, and boundaries
u0 = [0.0, 0.0];
p  = [1.0, 100.0];
lb = [-5.0, -5.0];
ub = [5.0, 5.0];

# 3. Create the OptimizationProblem
prob = OptimizationProblem(rosenbrock, u0, p, lb=lb, ub=ub);

sol = solve(prob, HMSSolver());

println("Best solution found: ", sol.u)
println("Objective value: ", sol.objective)
```

#### Obtaining original solution

One can access the original solver solution using original property.

``` @repl optim

sol.original
```

## Passing Configuration

> **Note**: Currently, HMS-specific parameters should be passed directly to the `HMSSolver` constructor. Passing HMS configuration arguments directly to the `solve` function is not supported at the moment.

To customize the behavior of the solver — such as adjusting the tree structure config, population sizes, the sprout condition and others — instantiate `HMSSolver` with the desired keywords.

### Example: Customized HMS Setup

In this example, we set a 3-level tree configuration with specific population sizes, sigma, global stop condition and `NormalPopulationCreator` population strategy to solve `eggholder` function.

``` @repl optim

function eggholder(x::Vector{Float64}, p)
    return -(x[2] + 47) * sin(sqrt(abs(x[2] + x[1]/2 + 47))) - x[1] * sin(sqrt(abs(x[1] - (x[2] + 47))))
end;

eggholder_problem = OptimizationProblem(
    eggholder, 
    [450.0, 400.0], 
    p, 
    lb=[-512.0, -512.0], 
    ub=[512.0, 512.0]
);

sigma = [[100.0, 100.0], [30.0, 30.0], [5.0, 5.0]];
level_config = [
    TreeLevelConfig(EvolutionaryGAMetaepoch),
    TreeLevelConfig(EvolutionaryGAMetaepoch),
    TreeLevelConfig(EvolutionaryCMAESMetaepoch),
];

# Configure the solver
solver = HMSSolver(
    level_config=level_config,
    population_sizes=[20, 10, 5],
    sigma=sigma,
    create_population=NormalPopulationCreator(),
    gsc=ProblemEvaluationLimitReached(5000) 
);

sol = solve(eggholder_problem, solver);
sol.original
```

``` @docs
HMSSolver
```
