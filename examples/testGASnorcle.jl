include("../src/MT1DGA.jl")
using MT1DGA

function testGASnorcle(size=100, gen=2000)
    inDataFile = "SNO-96-106.data"
    outModelFile = replace(inDataFile, ".data", ".genetic.model")
    myData = readdlm(inDataFile)

    N = 4
    zBounds = [LayerBC(0,0)
               fill(LayerBC(1,100_000), (N-1,1))]
    rBounds = fill(LayerBC(1,50_000), (N,1))

    pop = Population(size, myData, zBounds, rBounds)
    evolve!(pop, gen)

    writedlm(outModelFile, pop.cs[1].model)
    return pop
end
