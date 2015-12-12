include("../src/MT1DGeneticInversion.jl")
using MT1DGeneticInversion

function testGASnorcle(size=100, gen=2000)
    inDataFile = "SNO-96-106.data"
    outModelFile = replace(inDataFile, ".data", ".genetic.model")
    data = readdlm(inDataFile)

    N = 4
    zBounds = [LayerBC(0,0)
               fill(LayerBC(1,100_000), (N-1,1))]
    rBounds = fill(LayerBC(1,50_000), (N,1))

    I = Inversion(data, size, zBounds, rBounds)
    evolve!(I, gen)

    println("Best model:")
    println(I.pop[1].model)
    return I
end
