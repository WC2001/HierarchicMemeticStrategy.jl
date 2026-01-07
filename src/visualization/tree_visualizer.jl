struct HMSTreeVisualizer
    root::Deme
    demes::Vector{Deme}
    best::FitnessResult
end

function printTree(v::HMSTreeVisualizer; prefix::String = "", islast::Bool = true)
    root = v.root
    demes = v.demes
    best_sol = v.best.solution
    best_fit = v.best.fitness

    function adaptive_digits(x)
        absx = abs(x)
        if absx >= 0.01
            return 2
        elseif absx > 0
            return min(8, max(2, ceil(Int, -log10(absx)) + 2))
        else
            return 2
        end
    end

    function digits_for_vector(vec)
        if isempty(vec)
            return 2
        end
        maximum(adaptive_digits.(vec))
    end

    function fmt(x, digits)
        rounded = round(x, digits=digits)
        s = string(rounded)
        if !occursin('.', s)
            s *= "." * "0"^digits
        else
            parts = split(s, '.')
            decimal_len = length(parts[2])
            if decimal_len < digits
                s *= "0"^(digits - decimal_len)
            end
        end
        return s
    end

    function format_vector(vec::Vector{Float64})
        d = digits_for_vector(vec)
        return "(" * join(fmt.(vec, d), ", ") * ")"
    end

    function is_approx_equal(v1::Vector{Float64}, v2::Vector{Float64}; atol=1e-8)
        length(v1) == length(v2) || return false
        all(abs.(v1 .- v2) .<= atol)
    end

    function _printTree(node::Deme, prefix::String, islast::Bool)
        branch = node.parent_id === nothing ? "" : islast ? "└─ " : "├─ "
        active_str = node.is_active ? ", active" : ""
    
        if node.evaluations_count == 0  
            sprout_str = node.sprout === nothing ? "" : " | sprout:" * format_vector(node.sprout)
            println(prefix * branch, "(new deme)", active_str, sprout_str)
        else
            d_sol = digits_for_vector(node.best_solution)
            sol_str = "(" * join(fmt.(node.best_solution, d_sol), ", ") * ")"
    
            d_fit = adaptive_digits(node.best_fitness)
            fitness_str = fmt(node.best_fitness, d_fit)
    
            sprout_str = node.sprout === nothing ? "" : " | sprout:" * format_vector(node.sprout)
    
            star = is_approx_equal(node.best_solution, best_sol) ? " ★" : ""
    
            println(prefix * branch,
                "f", sol_str, " = ", fitness_str,
                star,
                ", evaluations: ", node.evaluations_count,
                active_str,
                sprout_str
            )
        end
    
        children = filter(d -> d.parent_id == node.id, demes)
        n = length(children)
        newprefix = prefix * (islast ? "   " : "│  ")
    
        for (i, child) in enumerate(children)
            _printTree(child, newprefix, i == n)
        end
    end
    

    _printTree(root, prefix, islast)
end


