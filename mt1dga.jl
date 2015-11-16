include("./mt1d.jl")
module MT1DGA
using GeneticAlgorithms
using MT1D
using Debug

global nLayers  = 3
global maxDepth = 100_000
global resRange = 1:100
global data     = zeros(1,5)
global mutationRate = 0.01

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
    global data
    fs = data[:,1].^-1
    ρa, Φ = mt1d(ent.model, fs)

    ρRMS = rms(data[:,2], ρa, data[:,4])
    ΦRMS = rms(data[:,3], Φ, data[:,5])
    totalRMS = sqrt(ρRMS^2 + ΦRMS^2)
end

# Scores closer to zero are better so override this.
function isless(lhs::MTModel, rhs::MTModel)
    abs(lhs.fitness) > abs(rhs.fitness)
end

function group_entities(pop::Array{MTModel})
    # This'll select the best models to breed
    # TODO: Add stopping conditions here
end

# TODO: Add type signature
@debug function crossover(parents)
    # Breeds two elite parents. The first layer of the child is
    # selected randomly from the first layer of the two models. Then,
    # a random pivot point is chosen. The layers above the pivot point
    # are taken from one random parent while the layers below are
    # taken from another (could be the same one). If the pivot point
    # is the last layer, all the layers beneath the first layer come
    # from the randomly selected above pivot parent.
    global nLayers
    childModel = zeros(nLayers, 2)

    topParent = rand(1:2)
    upperPivotParent = rand(1:2)
    lowerPivotParent = rand(1:2)

    pivotPoint = rand(2:nLayers)
    childModel[1,:] = parents[topParent].model[1,:]
    if pivotPoint == nLayers
        childModel[2:end,:] = parents[upperPivotParent].model[2:end,:]
    else
        childModel[2:pivotPoint,:] = parents[upperPivotParent].model[2:pivotPoint,:]
        childModel[pivotPoint+1:end,:] = parents[lowerPivotParent].model[pivotPoint+1:end,:]
    end

    MTModel(childModel)
end

function mutate(ent::MTModel)
    # Idea to mutate: Select a random layer to mutate the depth with,
    # then select a random layer to mutate the resistivity with.
    # How to mutate? -- For now, we'll just generate a new random number.
    # TODO: Explore alternate mutations
    global mutationRate, maxDepth, resRange;
    if rand(Float64) > mutationRate
        return
    end

    modelLayers = size(ent.model)[1]
    mutatedResLayer = rand(2:modelLayers)
    mutatedDepthLayer = rand(2:modelLayers)

    mutatedRes = rand(resRange)
    mutatedDepth = rand(1:maxDepth)

    # TODO: Check we don't exceed maxDepth and resRange
    ent.model[mutatedDepthLayer, 1] += mutatedDepth
    ent.model[mutatedResLayer, 2] += mutatedRes
end

end
