using PlotlyJS
using Colors: distinguishable_colors, hex
using DefaultApplication
using Dates

struct HMSResultVisualizer
    summaries::Vector{MetaepochSummary}
end

function plotDemeHistory(
    visualizer::HMSResultVisualizer, 
    prob::HMSOptimizationProblem,
    deme_index::Int, 
    x_index::Int, 
    y_index::Int
)
    summaries = visualizer.summaries
    n_summaries = length(summaries)

    if length(summaries[end].demes) < deme_index
        error("deme_index ($deme_index) is out of bounds; only $(length(summaries[end].demes)) demes created.")
    end

    bounds = problem_bounds(prob)
    genome_length = length(bounds.lower)

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

    x_range = [bounds.lower[x_index], bounds.upper[x_index]]
    y_range = use_zero_y ? [-1.0, 1.0] : [bounds.lower[y_index], bounds.upper[y_index]]

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
        xaxis = attr(title = "Dimension $x_index", range = x_range, autorange = false),
        yaxis = attr(title = "Dimension $y_index", range = y_range, autorange = false)
    )
    
    plot(traces, layout)
end

function plotPopulations(
    visualizer::HMSResultVisualizer,
    prob::HMSOptimizationProblem, 
    x_index::Int, 
    y_index::Int,
    filename::String = ""
)
    summaries = visualizer.summaries
    n_metaepochs = length(summaries)
    bounds = problem_bounds(prob)
    genome_length = length(bounds.lower)

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

    x_range = [bounds.lower[x_index], bounds.upper[x_index]]
    y_range = use_zero_y ? [-1.0, 1.0] : [bounds.lower[y_index], bounds.upper[y_index]]

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
                y_vals = use_zero_y ? zeros(length(genomes)) : [g[y_index] for g in genomes]
                
                marker_symbols = fill("circle", length(genomes))
                if epoch_idx == n_metaepochs && n_metaepochs > 1
                    prev_demes = summaries[n_metaepochs - 1].demes
                    if deme_idx <= length(prev_demes)
                        prev_count = length(prev_demes[deme_idx].population.genomes)
                        if length(genomes) == prev_count + 1
                            marker_symbols[end] = "x"
                        end
                    end
                end

                push!(deme_scatters,
                      scatter(
                          x = x_vals,
                          y = y_vals,
                          mode = "markers",
                          name = "Deme $deme_idx",
                          marker = attr(color = colors[deme_idx], opacity=1.0, symbol = marker_symbols)
                      ))
            end
        end

        epoch_traces[epoch_idx] = deme_scatters
    end

    initial_traces = epoch_traces[1]

    for traces in epoch_traces
        for (j, trace) in enumerate(traces)
            trace[:uid] = "trace-$j"
        end
    end

    frames = [
        frame(
            data = epoch_traces[i],
            name = "metaepoch$i",
            traces = collect(0:length(epoch_traces[1])-1)
        )
        for i in 1:n_metaepochs
    ]

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

    updatemenus = [
        attr(
            type = "buttons",
            buttons = [
                attr(
                    label = "Play",
                    method = "animate",
                    args = [nothing, attr(
                        frame = attr(duration = 500, redraw = true),
                        fromcurrent = true,
                        transition = attr(duration = 300, easing = "quadratic-in-out")
                    )]
                ),
                attr(
                    label = "Pause",
                    method = "animate",
                    args = [[nothing], attr(
                        mode = "immediate",
                        frame = attr(duration = 0, redraw = true),
                        transition = attr(duration = 0)
                    )]
                )
            ],
            direction = "left",
            pad = attr(r = 10, t = 87),
            showactive = false,
            x = 0.5,
            xanchor = "center",
            y = -0.2,
            yanchor = "top"
        )
    ]

    layout = Layout(
        title = "Populations by Metaepoch",
        xaxis = attr(title = "Dimension $x_index", range = x_range, autorange = false),
        yaxis = attr(title = use_zero_y ? "(0)" : "Dimension $y_index", range = y_range, autorange = false),
        sliders = sliders,
        updatemenus = updatemenus
    )
    
    p = Plot(initial_traces, layout, frames)
    
    final_filename = if isempty(filename)
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
        "hms_plot_$(timestamp).html"
    else
        endswith(filename, ".html") ? filename : filename * ".html"
    end

    savefig(p, final_filename)
    DefaultApplication.open(final_filename)
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

