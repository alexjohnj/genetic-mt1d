"""
Description
===========

The `Model` type represents a solution in the genetic algorithm.

Fields
======

- `model::Matrix`: The underlying model in matrix form. Each row is a layer with
  depth to the layer's top in the first column and the layer's resistivity in
  the second.
- `N::Integer`: The number of layers in the model.
- `zCodeParams::Array{LayerBC}`: The depth constraints for each layer.
- `rCodeParams::Array{LayerBC}`: The resistivity constraints for each layer.
- `fitness::Real`: The fitness of the model in the genetic algorithm. A -ve
  fitness indicates the fitness needs to be recalculated.
"""
type Model
    model::Matrix
    N::Integer
    zCodeParams::Array{LayerBC}
    rCodeParams::Array{LayerBC}
    fitness::Real

    function Model(m::Matrix, z::Array{LayerBC}, r::Array{LayerBC})
        if (z[1].min != 0 || z[1].max != 0)
            error("Range of depths for first layer must be (0,0)!")
        end
        if size(m)[1] != length(z) || size(m)[1] != length(r)
            error("Dimensions of zParams or rParams does not match number of layers")
        end

        new(m, size(m)[1], z, r, -1)
    end
end

"""
Description
===========

Create an instance of the `Model` type with a randomly generated model within a
set of constraints.

Arguments
=========

- `N::Integer`: The number of layers to use in the model.
- `zParams::Array{LayerBC}`: The depth constraints for each layer.
- `rParams::Array{LayerBC}`: The resistivity constraints for each layer.

Returns
=======

- `::Model`: An instance of the `Model` type.

Notes
=====

- **Throws an error** if `length(zParams) != length(rParams)`.
- **Throws an error** if the first element of `zParams != LayerBC(0,0)`.
"""
function createRandomModel(N::Integer, zParams::Array{LayerBC}, rParams::Array{LayerBC})
    if zParams[1].min != 0 || zParams[1].max != 0
        error("Range of depths for first layer must be (0,0)!")
    end

    if length(zParams) != N || length(rParams) != N
        error("Dimensions of zParams or rParams does not match number of layers")
    end

    model = zeros(N,2)
    model[1,2] = rand(rParams[1].min:rParams[1].max)
    for n = 2:N
        layerZParams = zParams[n]
        layerRParams = rParams[n]
        model[n,1] = rand(layerZParams.min:layerZParams.max)
        model[n,2] = rand(layerRParams.min:layerRParams.max)
    end

    Model(model, zParams, rParams)
end

"""
Description
===========

Calculate the fitness of a `Model` instance by calculating the RMS misfit
between the model's forward response and some data.

Arguments
=========

- `c::Model`: `Model` Instance to find the fitness of. The instance's `fitness`
  field will be mutated.
- `data::Matrix`: An N×5 matrix containing the data to find the misfit of. The
  column content should be:
    1. Period (s)
    2. Apparent resistivity (Ωm)
    3. Phase (°)
    4. Apparent resistivity uncertainty
    5. Phase uncertainty

Side Effects
============

- Updates the `fitness` field of the instance `c::Model` with the newly
  calculated misfit.

Returns
=======

- `nothing::Void`
"""
function calculateFitness!(c::Model, data::Matrix)
    if c.fitness >= 0
        # Fitness is already calculated
        return
    end

    fs = data[:,1].^-1
    ρ, Φ = calculateResponse(c.model, fs)

    RMSρ = calculateRMSMisfit(data[:,2], ρ, data[:,4])
    RMSΦ = calculateRMSMisfit(data[:,3], Φ, data[:,5])

    c.fitness = sqrt(RMSρ^2 + RMSΦ^2)

    # Are the model depths sorted from lowest to highest? Apply a
    # penalty if they aren't
    if !issorted(c.model[:,1])
        c.fitness += 50
    end

    # Do any of the layers in the current model exceed the search
    # boundaries? Apply a penalty proportional to how much they exceed
    # them if they do.
    zPenaltyFactor = 2
    rPenaltyFactor = 2
    for N in 1:c.N
        zCodeParams = c.zCodeParams[N]
        rCodeParams = c.rCodeParams[N]
        layerModel = c.model[N,:]

        if layerModel[1] < zCodeParams.min
            c.fitness += abs(zCodeParams.min - layerModel[1]) * zPenaltyFactor
        elseif layerModel[1] > zCodeParams.max
            c.fitness += abs(layerModel[1] - zCodeParams.max) * zPenaltyFactor
        end
        if layerModel[2] < rCodeParams.min
            c.fitness += abs(rCodeParams.min - layerModel[2]) * rPenaltyFactor
        elseif layerModel[2] > rCodeParams.max
            c.fitness += abs(layerModel[2] - rCodeParams.max) * rPenaltyFactor
        end
    end

    return
end


"""
Description
===========

Crossover two parent `Model` instances to produce two children using SBX.

Arguments
=========

- `a::Model`: The first parent.
- `b::Model`: The second parent.

Returns
=======

- `Tuple{Model,Model}`: The two offsprings of the parents.

Notes
=====

- The offsprings' fitness fields will not be initialised.
"""
function crossover(a::Model, b::Model, η::Integer)
    if η < 1
        error("Distribution parameter (η) must be a non-negative integer")
    end

    cA = deepcopy(a)
    cB = deepcopy(b)

    # Void the RMS calculations for the children
    cA.fitness = -1
    cB.fitness = -1

    for n = 1:a.N
        layerZParams = a.zCodeParams[n]
        layerRParams = a.rCodeParams[n]
        # Crossover depths
        if a.model[n,1] < b.model[n,1]
            βa, βb = boundedSBX(a.model[n,1], b.model[n,1], η, layerZParams.min, layerZParams.max)
            cA.model[n,1] = 0.5*(1+βa)*a.model[n,1] + 0.5 * (1-βa) * b.model[n,1]
            cB.model[n,1] = 0.5*(1-βb)*a.model[n,1] + 0.5 * (1+βb) * b.model[n,1]
        end

        # Crossover resistivities
        if a.model[n,2] < b.model[n,2]
            βa, βb = boundedSBX(a.model[n,2], b.model[n,2], η, layerRParams.min, layerRParams.max)
            cA.model[n,2] = 0.5*(1+βa)*a.model[n,2] + 0.5 * (1-βa) * b.model[n,2]
            cB.model[n,2] = 0.5*(1-βb)*a.model[n,2] + 0.5 * (1+βb) * b.model[n,2]
        end
    end

    # Force first layer to have zero depth
    cA.model[1,1] = 0
    cB.model[1,1] = 0

    cA.model = sortrows(cA.model)
    cB.model = sortrows(cB.model)

    return (cA, cB)
end

"""
Description
===========

Calculate the β values for the bounded SBX of two values.

Arguments
=========

- `p1::Real`: First parent value to crossover.
- `p2::Real`: Second parent value to crossover.
- `η::Integer`: Distribution parameter for the SBX.
- `a::Real`: Lower bound of search space.
- `b::Real`: Upper bound of search space.

Returns
=======

- `Tuple{Real,Real}`: The calculated values of β(u,a) and β(u,b).

Notes
=====

- **Throws an error** if `a > b` or if `p1 > p2`.
- Implementation is based on report by Deb, K. and Jain, H. (#2011017)
"""
function boundedSBX(p1::Real, p2::Real, η::Integer, a::Real, b::Real)
    if a > b
        error("boundedSBX: lower bound a was greater than b")
    end

    if !(p1 < p2)
        error("boundedSBX: Parent 1 was greater than parent 2")
    end
    βa = 1 + (p1 - a) / (p2 - p1)
    βb = 1 + (b - p2) / (p2 - p1)
    γa = 1 / (2 * βa^(η+1))
    γb = 1 / (2 * βb^(η+1))

    βua, βub = 0, 0

    ua = rand()
    if ua <= 0.5/(1-γa)
        βua = (2ua*(1-γa))^(1/(η+1))
    else
        βua = (1 / (2*(1-ua*(1-γa))))^(1/(η+1))
    end

    ub = rand()
    if ub <= 0.5/(1-γb)
        βub = (2ub*(1-γb))^(1/(η+1))
    else
        βub = (1 / (2*(1-ub*(1-γb))))^(1/(η+1))
    end

    return (βua, βub)
end

"""
Description
===========

Mutate a `Model` instance. For each depth and resistivity value in the model,
pick a random number. If it is less than the probability of mutation, pick a new
value within the constraints of the model.

Arguments
=========

- `c::Model`: Model to mutate.
- `Pm::Real`: Probability of mutation. Should be in the range [0,1].

Side Effects
============

- If a mutation occurs, `c`'s `model` field will be altered.

Returns
=======

- `nothing::Void`
"""
function mutate!(c::Model, Pm::Real)
    for n in 1:c.N
        # Mutate depth
        if rand() < Pm
            # Do nothing to the depth if this is the top layer
            if n != 1
                zParams = c.zCodeParams[n]
                c.model[n,1] = rand(zParams.min:zParams.max)
                sortrows(c.model)
                c.fitness = -1
            end
        end

        # Mutate resistivity
        if rand() < Pm
            rParams = c.rCodeParams[n]
            c.model[n,2] = rand(rParams.min:rParams.max)
            sortrows(c.model)
            c.fitness = -1
        end
    end

    return
end
