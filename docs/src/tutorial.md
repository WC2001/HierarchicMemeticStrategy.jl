```@meta
CurrentModule = hms_jl
```

```@setup result
using hms_jl
```

## Optimization

To demonstrate the package's capabilities, let's apply it to the
[Rosenbrock function](https://www.sfu.ca/~ssurjano/rosen.html),
a well-known benchmark in numerical optimization

### Defining a Problem

To optimize a function using **hms_jl**, you must encapsulate your objective function and its constraints within a FunctionProblem object.

First define the Rosenbrock function:

```@repl rosen
using hms_jl
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2;
```

Next, wrap it in a `FunctionProblem` object, specifying the bounds for each variable and whether to maximize or minimize:

```@repl rosen
lower_bounds = [-5.0, -5.0];
upper_bounds = [5.0, 5.0];

problem = FunctionProblem(
    fitness_function = rosenbrock,
    lower = lower_bounds,
    upper = upper_bounds,
    maximize = false
);
```
```@docs
FunctionProblem
```

### Performing Optimization

Once Function Problem is defined, we can perform optimization using hms algorithm.

```@repl rosen
res = hms(optimization_problem=problem)
```

### Configuration

The example above uses default configuration. For more configuration details 
check [HMS Configuration](@ref) section.

### Obtaining Results

Once the optimization process is complete, the solver returns an HMSResult object. This object contains the full history of the run, the best found solution, and tools for visualization.


#### List of functions

`HMSResult` interface for optimization result.

```@docs
HMSResult
summary(::HMSResult)
solution(::HMSResult)
best_fitness(::HMSResult)
iterations(::HMSResult)
f_calls(::HMSResult)
metaepoch_data(::HMSResult)
plotPopulations(::HMSResult, ::Int, ::Int)
plotBestFitness(::HMSResult)
plotDeme(::HMSResult, ::Int, ::Int, ::Int)
```



