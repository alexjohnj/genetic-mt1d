push!(LOAD_PATH, "../src")
using MT1DGeneticInversion
using MT1D
using Gadfly

"""
Forward models the apparent resistivity and phase of a 1D model and plots the
results using Gadfly. The plot is saved to `forward-modelled-response.html`.

## Arguments:
- `model` -- A matrix with depths in the first column and resistivities in the
  second. Ignore to use a default model.
"""
function testForwardModel(model=[])
    fs = logspace(-3, 3, 1024)
    if isempty(model)
        model = [0. 10;
                 1_000 100]
    end

    println("Forward modelling...")
    @time ρ, Φ = calc_response(model, fs)
    println("Done!")
    pRes = plot(x=fs.^-1,
                y=ρ,
                Geom.line,
                Guide.xlabel("Period (S)"),
                Guide.ylabel("ρa (Ωm)"),
                Guide.title("Apparent Resistivity"),
                Scale.x_log10,
                Scale.y_log10)
    pPhase = plot(x=fs.^-1,
                  y=Φ,
                  Geom.line,
                  Guide.xlabel("Period (S)"),
                  Guide.ylabel("Φ (°)"),
                  Guide.title("Phase"),
                  Scale.x_log10,
                  Scale.y_continuous(minvalue=0,maxvalue=90))
    pModel = plot(y=[model[:,2]; model[end,2]],
                  x=[model[:,1] / 1000; model[end,1] / 1000 + 0.5],
                  Stat.step(direction=:hv),
                  Geom.line,
                  Guide.ylabel("Resistivity (Ωm)"),
                  Guide.xlabel("Depth (km)"),
                  Guide.title("Model"),
                  Scale.y_log10)

    draw(SVG("forward-modelled-response.svg", 800px, 800px), hstack(vstack(pRes, pPhase), pModel))
    return
end

"""
Invert the data in the file ./SNO-96-106.data using the default inversion
parameters. This is data collected over the northern Canadian shield as part of
Lithoprobes's SNORCLE project. The data is the average of the TE and TM
measurements which were similar and so the structure is one-dimensional. The
inversion uses 4 layers and should produce a result with a drop in resistivity
at ~35 km, corresponding to the Moho discontinuity. See "The Electric Moho,
Jones & Ferguson (2001)" for a detailed analysis of this data set.

Arguments
=========

- `popsize::Integer`: Population size for the inversion
- `ngen::Integer`: Number of generations for inversion

Returns
=======

Nothing
"""
function testGASnorcle(popsize=150, ngen=3000)
    inDataFile = "./SNO-96-106.data"
    outDataFile = "./SNO-96-106.genetic.model"
    data = readdlm(inDataFile)

    nLayers = 4
    zBounds = [LayerBC(0, 0);
               fill(LayerBC(0, 100_000), (nLayers-1, 1))]
    rBounds = fill(LayerBC(0, 50_000), (nLayers, 1))

    I = Inversion(data, popsize, zBounds, rBounds)
    evolve!(I, ngen, verbose=true)
    print("Best model:\n$(I.pop[1].model)\n")
    p = plotResponseDataCompare(data, I.pop[1].model)
    draw(SVG("inversion-results.svg", 800px, 800px), p)
    return
end

function plotResponseDataCompare(data::Matrix, model::Matrix)
    # Calculate responses
    ρ, Φ = calc_response(model, data[:,1].^-1)

    resMin = data[:,2] - data[:,4]
    resMax = data[:,2] + data[:,4]
    dataResLayer = layer(x=data[:,1],
                         y=data[:,2],
                         ymin=resMin,
                         ymax=resMax,
                         Geom.point,
                         Geom.errorbar)
    phaseMin = data[:,3] - data[:,5]
    phaseMax = data[:,3] + data[:,5]
    dataPhaseLayer = layer(x=data[:,1],
                           y=data[:,3],
                           ymin=phaseMin,
                           ymax=phaseMax,
                           Geom.point,
                           Geom.errorbar)


    modelTheme = Theme(default_color=colorant"red")
    modelResLayer = layer(x=data[:,1],
                          y=ρ,
                          Geom.line,
                          modelTheme)
    modelPhaseLayer = layer(x=data[:,1],
                            y=Φ,
                            Geom.line,
                            modelTheme)

    resPlot = plot(dataResLayer,
                   modelResLayer,
                   Guide.xlabel("Period (S)"),
                   Guide.ylabel("ρa (Ωm)"),
                   Guide.title("Apparent Resistivity"),
                   Scale.x_log10,
                   Scale.y_log10)
    phasePlot = plot(dataPhaseLayer,
                     modelPhaseLayer,
                     Guide.xlabel("Period (S)"),
                     Guide.ylabel("Φ (°)"),
                     Guide.title("Phase"),
                     Scale.x_log10,
                     Scale.y_continuous(minvalue=0, maxvalue=90))

    modelPlot = plot(y=[model[:,2]; model[end,2]],
                     x=[model[:,1] / 1000; model[end,1] / 1000 + 1],
                     Stat.step(direction=:hv),
                     Geom.line,
                     Guide.ylabel("Resistivity (Ωm)"),
                     Guide.xlabel("Depth (km)"),
                     Guide.title("Model"),
                     Scale.y_log10,
                     modelTheme,
                     Guide.manual_color_key("Legend", ["Data", "Model"], ["deepskyblue", "red"]))

    return hstack(vstack(resPlot, phasePlot), modelPlot)
end
