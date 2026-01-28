# HMS Configuration

HMS algorithm allows for customizing numerous configuration elements.

```@docs
hms
```

## HMS Tree Configuration

To define the evolutionary algorithms at different levels of the tree, use `TreeLevelConfig`.

```@docs
TreeLevelConfig
DEFAULT_LEVEL_CONFIG
```

## HMS Evolutionary algorithms used in metaepochs

Here there are details for what evolutionary algorithms hms uses.

### MetaepochRunner Interface

All metaepochs must implement simple MetaepochRunner interface.

```@docs
MetaepochRunner
MetaepochResult
```

At present, `HierarchicMemeticStrategy.jl` provides 3 different metaepochs: 
- [`EvolutionaryGAMetaepoch`](@ref)
- [`EvolutionaryCMAESMetaepoch`](@ref)
- [`EvolutionaryDEMetaepoch`](@ref)

These metaepochs are built upon the robust implementations of [Genetic Algorithms (GA)](https://wildart.github.io/Evolutionary.jl/stable/ga/), [CMA-ES](https://wildart.github.io/Evolutionary.jl/stable/cmaes/), and [Differential Evolution (DE)](https://wildart.github.io/Evolutionary.jl/stable/de/) provided by the [Evolutionary.jl](https://github.com/wildart/Evolutionary.jl) package.

### EvolutionaryGAMetaepoch

```@docs
EvolutionaryGAMetaepoch
run_metaepoch(::EvolutionaryGAMetaepoch, ::Function, ::Bounds, ::Vector{Vector{Float64}}, ::Bool,
::Dict{String, Any})
```

### EvolutionaryCMAESMetaepoch

```@docs
run_metaepoch(::EvolutionaryCMAESMetaepoch, ::Function, ::Bounds, ::Vector{Vector{Float64}}, ::Bool,
::Dict{String, Any})
```

### EvolutionaryDEMetaepoch

```@docs
run_metaepoch(::EvolutionaryDEMetaepoch, ::Function, ::Bounds, ::Vector{Vector{Float64}}, ::Bool,
::Dict{String, Any})
```

## Stop conditions

Stop conditions are configurable hms parameters. They are used to control the process.

### Global stop condition

The Global Stop Condition (GSC) defines the termination criteria for the entire HMS algorithm, such as reaching a maximum number of function evaluations. It ensures the overall search process concludes once the primary optimization goals or resource limits are met.

#### Available Global Stop Conditions
The following conditions can be passed to the hms solver to control the total execution:

```@docs
ProblemEvaluationLimitReached
GlobalMetaepochLimitReached
```

#### Default GSC

The default GSC is set to ProblemEvaluationLimitReached with 10000 fitness function evaluation limit.

```julia
GSC = ProblemEvaluationLimitReached(10000)
```

### Local stop condition

The Local Stop Condition (LSC) determines when an individual deme should stop its local search. This allows the HMS tree to dynamically reallocate resources by pruning stagnant demes and sprouting new ones in more promising areas.

#### Available Local Stop Conditions
The following conditions can be passed to the hms solver to control the demes:

```@docs
LocalProblemEvaluationLimitReached
AllChildrenStopped
MetaepochWithoutBestFitnessImprovement
```
#### Default LSC

The default LSC is set to MetaepochWithoutBestFitnessImprovement with 10 metaepochs 
without improvement limit.

```julia
DEFAULT_LSC = MetaepochWithoutBestFitnessImprovement(10)
```
### Sprout condition

Configurable strategy parameter determining when a deme should sprout a child deme.

#### Sprout condition interface

All sprout conditions inherit from `AbstractSproutCondition`.

```@docs
HierarchicMemeticStrategy.AbstractSproutCondition
```

#### Sprout condition based on metrics

The primary way to maintain diversity in HMS is to enforce a minimum distance between demes at the same level. This prevents multiple demes from converging to the same local optimum.

```@docs
HierarchicMemeticStrategy.MetricSproutCondition
```

#### Default Sprout Condition

The default condition uses euclidean distance with default max distances calculated based on `sigma` parameter.

```@docs
HierarchicMemeticStrategy.DefaultSproutCondition
```

```julia
function sprout_default_euclidean_distances(sigma::Vector{Vector{Float64}})
    sprouting_condition_distance_ratio = 0.6
    return [sum(s .* sprouting_condition_distance_ratio) for s in sigma]
end
```

### Sigma

Parameter for evolutionary algorithms and Gaussian mutation. 
It is specified for each dimension at each level of the hms tree.

#### Default sigma

```@docs
HierarchicMemeticStrategy.default_sigma
```
### Population sizes

Parameter specifies number of individuals in demes at each level.

```@docs
HierarchicMemeticStrategy.default_population_sizes
```

### Population Initialization Interface

The population initialization strategy determines how individuals are distributed within the search space, both at the start of the optimization and when new demes "sprout" from parents.

#### The Abstract Population Strategy

All population initializers must inherit from `AbstractPopulationCreator`. This ensures a consistent, type-safe contract across different initialization methods.

```@docs
HierarchicMemeticStrategy.AbstractPopulationCreator
HierarchicMemeticStrategy.reseed!
```

#### Built-in Strategies

Package provides 2 strategies. The default one is `DefaultPopulationCreator`.

```@docs
HierarchicMemeticStrategy.HierarchicMemeticStrategy.DefaultPopulationCreator
HierarchicMemeticStrategy.NormalPopulationCreator
```

### Local optimizer

Parameter sets local optimizer to perform search from best individuals at leaf demes 
after global stop condition is satisfied. 

```@docs
LocalOptimizer
```

#### Default local optimizer

Default local optimizer uses **[L-BFGS](https://julianlsolvers.github.io/Optim.jl/stable/#algo/lbfgs/)** implementation from the **[Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)** package.

