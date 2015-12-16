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
- `depth_lims::Array{LayerBC}`: The depth constraints for each layer.
- `res_lims::Array{LayerBC}`: The resistivity constraints for each layer.
- `fitness::Real`: The fitness of the model in the genetic algorithm. A -ve
  fitness indicates the fitness needs to be recalculated.
"""
type Model
    model::Matrix
    N::Integer
    depth_lims::Array{LayerBC}
    res_lims::Array{LayerBC}
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
- `depth_lims::Array{LayerBC}`: The depth constraints for each layer.
- `res_lims::Array{LayerBC}`: The resistivity constraints for each layer.

Returns
=======

- `::Model`: An instance of the `Model` type.

Notes
=====

- **Throws an error** if `length(depth_lims) != length(res_lims)`.
- **Throws an error** if the first element of `depth_lims != LayerBC(0,0)`.
"""
function rand_model(N::Integer, depth_lims::Array{LayerBC}, res_lims::Array{LayerBC})
    if depth_lims[1].min != 0 || depth_lims[1].max != 0
        error("Range of depths for first layer must be (0,0)!")
    end

    if length(depth_lims) != N || length(res_lims) != N
        error("Dimensions of depth_lims or res_lims does not match number of layers")
    end

    model = zeros(N,2)
    model[1,2] = rand(res_lims[1].min:res_lims[1].max)
    for n = 2:N
        layer_depth_lim = depth_lims[n]
        layer_res_lim = res_lims[n]
        model[n,1] = rand(layer_depth_lim.min:layer_depth_lim.max)
        model[n,2] = rand(layer_res_lim.min:layer_res_lim.max)
    end

    Model(model, depth_lims, res_lims)
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
function calcfitness!(c::Model, data::Matrix)
    if c.fitness >= 0
        # Fitness is already calculated
        return
    end

    fs = data[:,1].^-1
    ρ, Φ = calc_response(c.model, fs)

    RMSρ = calc_rms(data[:,2], ρ, data[:,4])
    RMSΦ = calc_rms(data[:,3], Φ, data[:,5])

    c.fitness = sqrt(RMSρ^2 + RMSΦ^2)

    # Are the model depths sorted from lowest to highest? Apply a
    # penalty if they aren't
    if !issorted(c.model[:,1])
        c.fitness += 50
    end

    # Do any of the layers in the current model exceed the search
    # boundaries? Apply a penalty proportional to how much they exceed
    # them if they do.
    depth_penalty = 2
    res_penalty = 2
    for N in 1:c.N
        depth_lims = c.depth_lims[N]
        res_lims = c.res_lims[N]
        layer = c.model[N,:]

        if layer[1] < depth_lims.min
            c.fitness += abs(depth_lims.min - layer[1]) * depth_penalty
        elseif layer[1] > depth_lims.max
            c.fitness += abs(layer[1] - depth_lims.max) * depth_penalty
        end
        if layer[2] < res_lims.min
            c.fitness += abs(res_lims.min - layer[2]) * res_penalty
        elseif layer[2] > res_lims.max
            c.fitness += abs(layer[2] - res_lims.max) * res_penalty
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

    ca = deepcopy(a)
    cb = deepcopy(b)

    # Void the RMS calculations for the children
    ca.fitness = -1
    cb.fitness = -1

    for n = 1:a.N
        layer_depth_lims = a.depth_lims[n]
        layer_res_lims = a.res_lims[n]
        # Crossover depths
        if a.model[n,1] < b.model[n,1]
            βa, βb = boundedsbx(a.model[n,1],
                                b.model[n,1],
                                η,
                                layer_depth_lims.min,
                                layer_depth_lims.max)
            ca.model[n,1] = 0.5*(1+βa)*a.model[n,1] + 0.5 * (1-βa) * b.model[n,1]
            cb.model[n,1] = 0.5*(1-βb)*a.model[n,1] + 0.5 * (1+βb) * b.model[n,1]
        end

        # Crossover resistivities
        if a.model[n,2] < b.model[n,2]
            βa, βb = boundedsbx(a.model[n,2],
                                b.model[n,2],
                                η,
                                layer_res_lims.min,
                                layer_res_lims.max)
            ca.model[n,2] = 0.5*(1+βa)*a.model[n,2] + 0.5 * (1-βa) * b.model[n,2]
            cb.model[n,2] = 0.5*(1-βb)*a.model[n,2] + 0.5 * (1+βb) * b.model[n,2]
        end
    end

    # Force first layer to have zero depth
    ca.model[1,1] = 0
    cb.model[1,1] = 0

    ca.model = sortrows(ca.model)
    cb.model = sortrows(cb.model)

    return (ca, cb)
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
function boundedsbx(p1::Real, p2::Real, η::Integer, a::Real, b::Real)
    if a > b
        error("boundedsbx: lower bound a was greater than b")
    end

    if !(p1 < p2)
        error("boundedsbx: Parent 1 was greater than parent 2")
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
                layer_depth_lims = c.depth_lims[n]
                c.model[n,1] = rand(layer_depth_lims.min:layer_depth_lims.max)
                sortrows(c.model)
                c.fitness = -1
            end
        end

        # Mutate resistivity
        if rand() < Pm
            layer_res_params = c.res_lims[n]
            c.model[n,2] = rand(layer_res_params.min:layer_res_params.max)
            sortrows(c.model)
            c.fitness = -1
        end
    end

    return
end
