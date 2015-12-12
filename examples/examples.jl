push!(LOAD_PATH, "../src")
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
    @time ρ, Φ = calculateResponse(model, fs)
    println("Done!")
    pRes = plot(x=fs.^-1,
                y=ρ,
                Geom.line,
                Guide.xlabel("Period (S)"),
                Guide.ylabel("Apparent Res (Ωm)"),
                Guide.title("Apparent Resistivity"),
                Scale.x_log10,
                Scale.y_log10)
    pPhase = plot(x=fs.^-1,
                  y=Φ,
                  Geom.line,
                  Guide.xlabel("Period (S)"),
                  Guide.ylabel("Phase (°)"),
                  Guide.title("Phase (°)"),
                  Scale.x_log10,
                  Scale.y_continuous(minvalue=0,maxvalue=90))

    draw(SVG("forward-modelled-response.svg", 600px, 600px), vstack(pRes, pPhase))
end
