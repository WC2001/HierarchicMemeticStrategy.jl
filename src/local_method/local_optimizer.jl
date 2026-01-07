"""
    LocalOptimizer

An abstract type for all local refinement strategies in the HMS algorithm.

# Description
`LocalOptimizer` is used to perform fine-grained refinement on the 
best individual of a leaf demes after global stop condition is satisfied.

# Interface
Any subtype must implement:
`optimize!(optimizer::LocalOptimizer, deme::Deme, cached_f::CustomCache)`
"""
abstract type LocalOptimizer end

function optimize!(optimizer::LocalOptimizer, deme::Deme, cached_f::CustomCache)
    throw(MethodError(optimize!, (optimizer, deme, cached_f)))
end