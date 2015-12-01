include("./mt1d.jl")
using MT1D

fs = logspace(-3,3,40)
modelFile = "./data/synthetic-6-layer.model"
outFile = replace(modelFile, ".model", ".data", 1)

model = readdlm(modelFile, Float64)
ρa, Φ = mt1d(model, fs)

#Add some noise to the signal
Φ = Φ + 0.1 * Φ .* randn(length(Φ), 1)
ρa = ρa + 0.1 * ρa .* randn(length(ρa), 1)

data = hcat(fs.^-1, ρa, Φ, ρa * 0.12, fill(3, (40,1)))
data = round(data, 6);
writedlm(outFile, data);
