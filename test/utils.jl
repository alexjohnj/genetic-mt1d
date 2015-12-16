module Utils
using Gadfly

export plot_response

"""
Plot a comparison between the expected and calculated apparent resistivity and
phase response as well as the accompanying model.

Arguments
=========

- `expected::Matrix`: Expected response, ρa in column 1 and Φ in column 2.
- `calculated::Matrix`: Calculated response, ""
- `model::Matrix`: Model used to find response. Depth in column 1, ρ in column 2.
- `fs::Vector`: Frequencies used in forward modelling calculations

Returns
=======

- `::Plot`: Stacked plots of resistivity, phase and the model.
"""
function plot_response(expected::Matrix, calculated::Matrix, model::Matrix, fs::Vector)
    expected_theme = Theme(default_color=colorant"orange")
    lres_expected = layer(x=fs,
                          y=expected[:,1],
                          Geom.point,
                          expected_theme)
    lres_calculated = layer(x=fs,
                            y=calculated[:,1],
                            Geom.line)

    res_plot = plot(lres_expected,
                    lres_calculated,
                    Guide.xlabel("Frequency (Hz)"),
                    Guide.ylabel("ρa (Ωm)"),
                    Scale.x_log10,
                    Scale.y_log10,
                    Coord.cartesian(xflip=true))

    lphase_expected = layer(x=fs,
                            y=expected[:,2],
                            Geom.point,
                            expected_theme)
    lphase_calculated = layer(x=fs,
                              y=calculated[:,2],
                              Geom.line)
    phase_plot = plot(lphase_expected,
                      lphase_calculated,
                      Guide.xlabel("Frequency (Hz)"),
                      Guide.ylabel("Degree (°)"),
                      Scale.x_log10,
                      Scale.y_continuous(minvalue=0, maxvalue=90),
                      Coord.cartesian(xflip=true))

    model_plot = plot(y=[model[:,2]; model[end,2]],
                      x=[model[:,1] / 1000; model[end,1] / 1000 + 0.5],
                      Stat.step(direction=:hv),
                      Geom.line,
                      Guide.ylabel("Resistivity (Ωm)"),
                      Guide.xlabel("Depth (km)"),
                      Guide.title("Model"),
                      Scale.y_log10,
                      Guide.manual_color_key("Legend", ["Expected";"Calculated"], ["orange";"deepskyblue"]))

    return hstack(vstack(res_plot, phase_plot), model_plot)
end

end
