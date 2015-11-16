include("./MT1DModel.jl")
module MT1DGA
using MT1DModel

export MT1DModelPopulation, sortPopulation!, matePopulation!

type MT1DModelPopulation
    size::Integer
    layers::Integer
    zMax::Integer
    resRange::Range
    data

    pop
    gen::Integer

    function MT1DModelPopulation(size, layers, zMax, resRange, data)
        P = []
        for i in 1:size
            push!(P, createRandomModel(layers, zMax, resRange, data))
        end
        new(size, layers, zMax, resRange, data, P, 1)
    end
end

function sortPopulation!(pop::Array)
    sort!(pop, lt=(a,b) -> a.RMS < b.RMS)
end

function cullPopulation!(MTPop::MT1DModelPopulation, elitismRate::Real)
    targetSize = floor(Integer, elitismRate * MTPop.size)
    splice!(MTPop.pop, targetSize+1:MTPop.size)
end

function matePopulation!(MTPop::MT1DModelPopulation)
    nParents = length(MTPop.pop)
    for i = 1:2:nParents-1
        child = crossover(MTPop.pop[i], MTPop.pop[i+1])
        push!(MTPop.pop, child)
    end
end

function fillPopulation!(MTPop::MT1DModelPopulation)
    popToMake = MTPop.size - length(MTPop.pop)
    for i in 1:popToMake
        push!(MTPop.pop, createRandomModel(MTPop.layers, MTPop.zMax, MTPop.resRange, MTPop.data))
    end
end

function advanceGeneration(MTPop::MT1DModelPopulation)
    sortPopulation!(MT1DModelPopulation.pop)
    cullPopulation!(MT1DModelPopulation.pop, 0.6)
    matePopulation!(MT1DModelPopulation.pop)
    fillPopulation!(MT1DModelPopulation.pop)
    sortPopulation!(MT1DModelPopulation.pop)
end

end
