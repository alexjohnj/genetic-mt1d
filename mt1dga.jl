include("./MT1DModel.jl")
module MT1DGA
using MT1DModel

export MT1DModelPopulation, evolve!

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
    sort!(pop, lt=(a,b) -> a.fitness < b.fitness)
end

function cullPopulation!(MTPop::MT1DModelPopulation, elitismRate::Real)
    targetSize = floor(Integer, elitismRate * MTPop.size)
    splice!(MTPop.pop, targetSize+1:MTPop.size)
end

function matePopulation!(MTPop::MT1DModelPopulation)
    nParents = length(MTPop.pop)
    for i = 1:2:nParents-1
        childA, childB = crossover(MTPop.pop[i], MTPop.pop[i+1])
        calculateFitness!(childA, MTPop.data)
        calculateFitness!(childB, MTPop.data)
        push!(MTPop.pop, childA, childB)
        if length(MTPop.pop) >= MTPop.size
            break
        end
    end
end

function fillPopulation!(MTPop::MT1DModelPopulation)
    popToMake = MTPop.size - length(MTPop.pop)
    for i in 1:popToMake
        push!(MTPop.pop, createRandomModel(MTPop.layers, MTPop.zMax, MTPop.resRange, MTPop.data))
    end
end

function mutatePopulation!(MTPop::MT1DModelPopulation, mutationRate::Real)
    for i in 1:MTPop.size
        mutate!(MTPop.pop[i], mutationRate)
    end
end

function advanceGeneration!(MTPop::MT1DModelPopulation)
    sortPopulation!(MTPop.pop)
    cullPopulation!(MTPop, 0.70)
    matePopulation!(MTPop)
    fillPopulation!(MTPop)
    mutatePopulation!(MTPop, 0.01)
    sortPopulation!(MTPop.pop)
    MTPop.gen += 1
end

function evolve!(MTPop::MT1DModelPopulation, maxgen=1000)
    while MTPop.gen < maxgen
        advanceGeneration!(MTPop)
        if MTPop.pop[1].fitness <= 1.5
            println("*"^80)
            @printf("Gen: %d\tFitness: %.2f\tSize: %d\n", MTPop.gen, MTPop.pop[1].fitness, length(MTPop.pop))
            break
        end
        if MTPop.gen % 50 == 0
            @printf("Gen: %d\tFitness: %.2f\tSize: %d\n", MTPop.gen, MTPop.pop[1].fitness, length(MTPop.pop))
        end
    end
end

end
