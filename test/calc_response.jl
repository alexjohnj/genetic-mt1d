"""
Tests the forward modelling function `calculate_response` against a range of
models. For model information, see "./test-models/README.md". A comparison
of the expected output and calculated output for each model is plotted and saved
into "./test-models/generated-plots".
"""
module TestCalcResponse
using Base.Test
using MT1D
using Utils
using Gadfly

fs = logspace(-4, 4, 100)
cd("./test-models") do
    model_files = filter(s -> ismatch(r"\.model$", s), readdir())

    map(model_files) do fname
        model = readdlm(fname)
        expected_response = readdlm("$fname.response")

        ρ, Φ = calc_response(model, fs)
        p = plot_response(expected_response, [ρ Φ], model, fs)
        draw(SVG("./generated-plots/$fname.svg", 1140px, 800px), p)

        @test_approx_eq expected_response[:,1] ρ
        @test_approx_eq expected_response[:,2] Φ
    end
end

end
