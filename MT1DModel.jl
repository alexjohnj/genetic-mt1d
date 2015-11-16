module MT1DModel

export createRandomModel

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

end
