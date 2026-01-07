using PlotlyJS
using Colors: distinguishable_colors, hex
using DefaultApplication

struct HMSResultVisualizer
    summaries::Vector{MetaepochSummary}
end

function plotDemeHistory(visualizer::HMSResultVisualizer, deme_index::Int, x_index::Int, y_index::Int)
    summaries = visualizer.summaries
    n_summaries = length(summaries)

    if length(summaries[end].demes) < deme_index
        error("deme_index ($deme_index) is out of bounds; only $(length(summaries[end].demes)) demes created.")
    end

    genome_length = length(summaries[1].demes[1].population.genomes[1])

    if x_index < 1 || x_index > genome_length
        error("x_index ($x_index) out of bounds for genome length $genome_length")
    end

    use_zero_y = false
    if y_index < 1 || y_index > genome_length
        if genome_length == 1
            use_zero_y = true
        else
            error("y_index ($y_index) out of bounds for genome length $genome_length")
        end
    end

    colors = distinguishable_colors(n_summaries)
    traces = GenericTrace[]
    prev_genome_count = nothing
    
    for (i, summary) in enumerate(summaries)
        if deme_index > length(summary.demes)
            continue
        end

        deme = summary.demes[deme_index]
        genomes = deme.population.genomes
        genome_count = length(genomes)
        
        x_vals = [genome[x_index] for genome in genomes]
        y_vals = use_zero_y ? zeros(length(genomes)) : [genome[y_index] for genome in genomes]

        symbols = fill("circle", genome_count)

        if !isnothing(prev_genome_count) && genome_count == prev_genome_count + 1
            symbols[end] = "x"
        end
        
        push!(traces, scatter(
            x = x_vals,
            y = y_vals,
            mode = "markers",
            name = "Metaepoch $i",
            marker = attr(color = colors[i], symbol = symbols)
        ))

        prev_genome_count = genome_count
    end
    
    layout = Layout(
        title = "Deme $deme_index",
        xaxis = attr(title = "x_$x_index"),
        yaxis = attr(title = "x_$y_index")
    )
    
    plot(traces, layout)
end


function savePopulationsPNGs(visualizer::HMSResultVisualizer, x_index::Int, y_index::Int; output_dir::String = "/Users/a.wiktor/Downloads/tmp_plots")
    summaries = visualizer.summaries
    n_metaepochs = length(summaries)
    genome_length = length(summaries[1].demes[1].population.genomes[1])

    if x_index < 1 || x_index > genome_length
        error("x_index ($x_index) out of bounds for genome length $genome_length")
    end

    use_zero_y = false
    if y_index < 1 || y_index > genome_length
        if genome_length == 1
            use_zero_y = true
        else
            error("y_index ($y_index) out of bounds for genome length $genome_length")
        end
    end

    n_demes = length(summaries[end].demes)
    colors = distinguishable_colors(n_demes)

    for epoch_idx in 1:n_metaepochs
        deme_scatters = GenericTrace[]

        for deme_idx in 1:n_demes
            deme = epoch_idx <= length(summaries) && deme_idx <= length(summaries[epoch_idx].demes) ?
                summaries[epoch_idx].demes[deme_idx] : nothing

            if deme === nothing
                push!(deme_scatters, scatter(x=[], y=[], mode="markers", marker=attr(opacity=0.0)))
            else
                genomes = deme.population.genomes
                x_vals = [genome[x_index] for genome in genomes]
                y_vals = use_zero_y ? zeros(length(genomes)) : [genome[y_index] for genome in genomes]

                push!(deme_scatters, scatter(
                    x = x_vals,
                    y = y_vals,
                    mode = "markers",
                    name = "Deme $deme_idx",
                    marker = attr(color = colors[deme_idx])
                ))
            end
        end

        layout = Layout(
            title = "Populations - Metaepoch $epoch_idx",
            xaxis = attr(title = "x_$x_index"),
            yaxis = attr(title = use_zero_y ? "(0)" : "x_$y_index"),
        )

        plt = plot(deme_scatters, layout)
        filename = joinpath(output_dir, "metaepoch_$(lpad(string(epoch_idx), 2, '0')).png")
        savefig(plt, filename)
    end
end



function plotPopulations(visualizer::HMSResultVisualizer, x_index::Int, y_index::Int)
    summaries = visualizer.summaries
    n_metaepochs = length(summaries)
    genome_length = length(summaries[1].demes[1].population.genomes[1])

    if x_index < 1 || x_index > genome_length
        error("x_index ($x_index) out of bounds for genome length $genome_length")
    end

    use_zero_y = false
    if y_index < 1 || y_index > genome_length
        if genome_length == 1
            use_zero_y = true
        else
            error("y_index ($y_index) out of bounds for genome length $genome_length")
        end
    end

    n_demes = length(summaries[end].demes)
    colors = distinguishable_colors(n_demes)

    
    epoch_traces = Vector{Vector{GenericTrace}}(undef, n_metaepochs)

    for epoch_idx in 1:n_metaepochs
        deme_scatters = GenericTrace[]

        for deme_idx in 1:n_demes
            deme = epoch_idx <= length(summaries) &&
                   deme_idx <= length(summaries[epoch_idx].demes) ?
                   summaries[epoch_idx].demes[deme_idx] : nothing

            if deme === nothing
                push!(deme_scatters, scatter(x=[], y=[], mode="markers",
                                             marker=attr(opacity=0.0)))
            else
                genomes = deme.population.genomes
                x_vals = [g[x_index] for g in genomes]
                y_vals = use_zero_y ? zeros(length(genomes)) :
                                      [g[y_index] for g in genomes]

                push!(deme_scatters,
                      scatter(
                          x = x_vals,
                          y = y_vals,
                          mode = "markers",
                          name = "Deme $deme_idx",
                          marker = attr(color = colors[deme_idx], opacity=1.0)
                      ))
            end
        end

        epoch_traces[epoch_idx] = deme_scatters
    end

    initial_traces = epoch_traces[1]

    # assign unique IDs to preserve color consistency
    for traces in epoch_traces
        for (j, trace) in enumerate(traces)
            trace[:uid] = "trace-$j"
        end
    end

    # --------------------
    # Create frames
    # --------------------
    frames = [
        frame(
            data = epoch_traces[i],
            name = "metaepoch$i",
            traces = collect(0:length(epoch_traces[1])-1)
        )
        for i in 1:n_metaepochs
    ]

    # --------------------
    # Slider (controls frames)
    # --------------------
    steps = [
        attr(
            method = "animate",
            args = [[frames[i].name],
                    attr(mode="immediate",
                         frame=attr(duration=0, redraw=true),
                         transition=attr(duration=0))],
            label = "Metaepoch $i"
        )
        for i in 1:n_metaepochs
    ]

    sliders = [attr(
        steps = steps,
        active = 0,
        currentvalue = attr(prefix = "Metaepoch: "),
        pad = attr(t = 50),
        x = 0, y = 0,
        len = 1.0
    )]

    # --------------------
    # Layout
    # --------------------
    layout = Layout(
        title = "Populations by Metaepoch",
        xaxis = attr(title = "x_$x_index"),
        yaxis = attr(title = use_zero_y ? "(0)" : "x_$y_index"),
        sliders = sliders
    )
    
    p = Plot(initial_traces, layout, frames)
    
    filename = "populations_x$(x_index)_y$(y_index).html"
    savefig(p, filename)
    DefaultApplication.open(filename)
end


function format_smart(value::Real)
    absval = abs(value)
    if absval < 1e-4
        return string(round(value, sigdigits=5))
    elseif absval < 0.01
        return strip_zeros(round(value, digits=4))
    elseif absval < 0.1
        return strip_zeros(round(value, digits=3))
    elseif absval < 1
        return strip_zeros(round(value, digits=2))
    elseif absval < 10
        return strip_zeros(round(value, digits=2))
    elseif absval < 100
        return strip_zeros(round(value, digits=1))
    else
        return strip_zeros(round(value))
    end
end

function strip_zeros(x::Real)
    s = string(x)
    s = replace(s, r"\.0+$" => "")              
    s = replace(s, r"(\.\d*?[1-9])0+$" => s"\1") 
    return s
end


function plotBestFitness(pv::HMSResultVisualizer)
    best_fitnesses = [summary.best_fitness for summary in pv.summaries]
    epochs = 1:length(best_fitnesses)

    trace = scatter(
        x = epochs,
        y = best_fitnesses,
        mode = "lines+markers",
        name = "Best fitness",
        line = attr(color = "blue"),
        marker = attr(color = "blue")
    )

    layout = Layout(
        title = "Best fitness over metaepochs",
        xaxis = attr(title = "Metaepoch"),
        yaxis = attr(title = "Best fitness"),
        hovermode = "closest"           
    )

    plot([trace], layout)
end

