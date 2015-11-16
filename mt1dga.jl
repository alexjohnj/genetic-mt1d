include("./MT1DModel.jl")
module MT1DGA
using MT1DModel

export MT1DModelPopulation, sortPopulation!

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

# function advanceGeneration(MTPop::MT1DModelPopulation)

# end

end
