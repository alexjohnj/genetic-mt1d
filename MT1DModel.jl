include("./mt1d.jl")
module MT1DModel
using MT1D

export createRandomModel, calculateRMS!

type MTModel
    model::Matrix
    N::Integer
    RMS

    function MTModel(model::Matrix)
        new(model, size(model)[1], nothing)
    end
end

"""
Creates a `MTModel` type with a model matrix of size Nx2 and no
calculated RMS.
# Arguments:
- `N`        -- The number of layers to use.
- `zMax`     -- The maximum depth to the top of the layers.
- `resRange` -- The range of resistivities to use in the model.

# Returns:
- An `MTModel` type.
"""
function createRandomModel(N::Integer, zMax::Integer, resRange::Range)
    model = zeros(N, 2)
    model = zeros(N, 2)
    model[1,2] = rand(resRange)
    model[2:end, 1] = rand(1:zMax, (N-1, 1))
    model[2:end, 2] = rand(resRange, (N-1, 1))

    MTModel(model)
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

    return totalRMS;
end

end
