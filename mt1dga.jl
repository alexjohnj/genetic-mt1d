include("./mt1d.jl")
module MT1DGA
using GeneticAlgorithms
using MT1D

global nLayers  = 3
global maxDepth = 100_000
global resRange = 1:100
global data     = zeros(1,5)

export nLayers, maxDepth, resRange, data
export setNLayers, setMaxDepth, setResRange, setData

function setNLayers(N::Integer)
    global nLayers
    nLayers = N
end

function setMaxDepth(z::Integer)
    global maxDepth
    maxDepth = z
end

function setResRange(rng::Range)
    global resRange
    resRange = rng
end

function setData(d::Matrix)
    global data
    data = d
end

type MTModel <: Entity
    model::Matrix
    fitness

    function MTModel(model::Matrix)
        new(model, nothing)
    end
end

function create_entity(idx::Integer)
    global nLayers, maxDepth, resRange
    model = zeros(nLayers, 2)
    model[1,2] = rand(resRange)
    model[2:end, 1] = rand(1:maxDepth, (nLayers-1, 1))
    model[2:end, 2] = rand(resRange, (nLayers-1, 1))

    MTModel(model)
end

function fitness(ent::MTModel)
    # Here we will calculate each entity's RMS misfit
end

# Scores closer to zero are better so override this.
function isless(lhs::MTModel, rhs::MTModel)
    abs(lhs.fitness) > abs(rhs.fitness)
end

function group_entities(pop::Array{MTModel})
    # This'll select the best models to breed
end

function crossover(ents::Array{MTModel})
    # This'll breed the best models somehow
end

function mutate(ent::MTModel)
    # This'll mutate a model, somehow
end

end
