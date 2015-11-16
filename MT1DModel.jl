include("./mt1d.jl")
module MT1DModel
using MT1D

export createRandomModel, calculateRMS!, mutate!

type MTModel
    model::Matrix
    N::Integer
    zMax::Integer
    resRange::Range
    RMS

    function MTModel(model::Matrix, zMax::Integer, resRange::Range)
        new(model, size(model)[1], zMax, resRange, nothing)
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
    model[2:end, 1] = rand(1:zMax, (N-1, 1))
    model[2:end, 2] = rand(resRange, (N-1, 1))

    MTModel(model, zMax, resRange)
end

function createRandomModel(N::Integer, zMax::Integer, resRange::Range, data)
    model = createRandomModel(N, zMax, resRange)
    calculateRMS!(model, data)
    return model
end

"""
Updates the `RMS` field of an MTModel type using the data
matrix. Matrix should be of the form, [T,ρ,Φ,σρ,σΦ].
## Arguments:
- `model` -- The `MTModel` to calculate the RMS of.
- `data`  -- The data matrix to use in the RMS calculation.

# Returns:
- The total RMS for the MTModel type.
"""
function calculateRMS!(model::MTModel, data)
    fs = data[:,1].^-1;
    ρa, Φ = mt1d(model.model, fs)

    ρRMS = rms(data[:,2], ρa, data[:,4])
    ΦRMS = rms(data[:,3], Φ, data[:,5])
    totalRMS = sqrt(ρRMS^2 + ΦRMS^2)

    model.RMS = totalRMS

    return totalRMS
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
    model.model[mutatedResLayer,2] += mutatedRes
end

end
