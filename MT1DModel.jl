include("./mt1d.jl")
module MT1DModel
using MT1D

export createRandomModel, calculateFitness!, mutate!, crossover

type MTModel
    model::Matrix
    N::Integer
    zMax::Integer
    resRange::Range
    fitness::Real

    function MTModel(model::Matrix, zMax::Integer, resRange::Range)
        new(model, size(model)[1], zMax, resRange, NaN)
    end
end

"""
Creates a `MTModel` type with a model matrix of size Nx2 and no
calculated RMS. If a data matrix is supplied, the model will also have
its RMS calculated.
# Arguments:
- `N`          -- The number of layers to use.
- `zMax`       -- The maximum depth to the top of the layers.
- `resRange`   -- The range of resistivities to use in the model.
- `data` (Opt) -- A data matrix from which to calculate an RMS misfit
  for the model.

# Returns:
- An `MTModel` type.
"""
function createRandomModel(N::Integer, zMax::Integer, resRange::Range)
    model = zeros(N, 2)
    model = zeros(N, 2)

    model[1,2] = rand(resRange)
    for n in 2:N
        depthRange = (n - 1) * zMax / N : n*zMax/N
        model[n, 1] = rand(depthRange)
    end
    model[2:end, 2] = rand(resRange, (N-1, 1))
    model = sortrows(model)

    MTModel(model, zMax, resRange)
end

function createRandomModel(N::Integer, zMax::Integer, resRange::Range, data)
    model = createRandomModel(N, zMax, resRange)
    calculateFitness!(model, data)
    return model
end

"""
Updates the `fitness` field of an MTModel type using the data
matrix. Matrix should be of the form, [T,ρ,Φ,σρ,σΦ]. Fitness is
defined as the quadratic norm of the model and data.
## Arguments:
- `model` -- The `MTModel` to calculate the RMS of.
- `data`  -- The data matrix to use in the RMS calculation.

# Returns:
- The current fitness of the model.
"""
function calculateFitness!(model::MTModel, data)
    fs = data[:,1].^-1;
    ρa, Φ = mt1d(model.model, fs)

    qNorm = norm(data[:,2] - ρa)^2 + norm(data[:,3] - Φ)^2
    model.fitness = qNorm
    return qNorm
end

"""
Selects two random layers to mutate their depth and resistivity.

## Arguments:
`model` -- The `MTModel` type to mutate.
`chance` -- The probability of mutation (between 0 and 1).
"""
function mutate!(model::MTModel, chance::Real)
    if rand() > chance
        return
    end

    mutatedResLayer = rand(2:model.N)
    mutatedDepthLayer = rand(2:model.N)

    mutatedRes = rand(model.resRange)
    mutatedDepth = rand(1:model.zMax)

    # TODO: Check we don't exceed maxDepth and resRange
    model.model[mutatedDepthLayer,1] += mutatedDepth
    while model.model[mutatedDepthLayer,1] > model.zMax
        model.model[mutatedDepthLayer,1] -= rand(1:model.zMax)
    end

    model.model[mutatedResLayer,2] += mutatedRes
    while model.model[mutatedResLayer,2] > maximum(model.resRange)
        model.model[mutatedResLayer,2] -= rand(model.resRange)
    end
    model.model = sortrows(model.model)
end

function crossover(modelA::MTModel, modelB::MTModel)
    # TODO: Audit this and maybe? rewrite
    # Select a random layer to swap between the models

    # swapLayer = rand(1:modelA.N)
    # childModelA = deepcopy(modelA)
    # childModelB = deepcopy(modelB)
    # childModelA.model[swapLayer,:], childModelB.model[swapLayer,:] = childModelB.model[swapLayer,:], childModelA.model[swapLayer,:]
    # childModelA.model = sortrows(childModelA.model)
    # childModelB.model = sortrows(childModelB.model)
    # return (childModelA, childModelB)

    parents = (modelA, modelB)

    childModelA = zeros(modelA.N, 2)
    childModelB = zeros(modelA.N, 2)

    childModelA[1,:] = modelB.model[1,:]
    childModelB[1,:] = modelA.model[1,:]

    pivotPoint = rand(2:modelA.N)
    if pivotPoint == modelA.N
        childModelA[2:end,:] = modelB.model[2:end,:]
        childModelB[2:end,:] = modelA.model[2:end,:]
    else
        childModelA[2:pivotPoint,:] = modelB.model[2:pivotPoint,:]
        childModelA[pivotPoint+1:end,:] = modelA.model[pivotPoint+1:end,:]
        childModelB[2:pivotPoint,:] = modelA.model[2:pivotPoint,:]
        childModelB[pivotPoint+1:end,:] = modelB.model[pivotPoint+1:end,:]
    end

    childModelA = sortrows(childModelA)
    childModelB = sortrows(childModelB)
    # Assume models have the same zMax and resRange
    (MTModel(childModelA, modelA.zMax, modelA.resRange), MTModel(childModelB, modelA.zMax, modelA.resRange))
end

end
