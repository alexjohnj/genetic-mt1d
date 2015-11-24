include("./mt1d.jl")
using MT1D

fs = 3 * logspace(-3,3,40)
modelFile = "./data/synthetic-4-layer.model"
outFile = replace(modelFile, ".model", ".data", 1)

model = readdlm(modelFile, Float64)
ρa, Φ = mt1d(model, fs)

data = hcat(fs.^-1, ρa, Φ, ρa * 0.12, fill(3, (40,1)))
writedlm(outFile, data);
