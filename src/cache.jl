using ConcurrentCollections

struct CustomCache{F}
    f::F
    cache::ConcurrentDict{Any, Any}
end

function CustomCache(f::F) where F
    CustomCache(f, ConcurrentDict{Any, Any}())
end

function (c::CustomCache)(x)
    if haskey(c.cache, x)
        return c.cache[x]
    else
        val = c.f(x)
        c.cache[x] = val
        return val
    end
end

function iscached(c::CustomCache, x)
    return haskey(c.cache, x)
end

