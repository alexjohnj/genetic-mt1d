include("./mt1d.jl")
using MT1D

type LayerBC
    min::Integer
    max::Integer
    res::Integer # Number of significant figures desired for parameter
end

type Chromosome
    model::Matrix
    N::Integer
    zCodeParams::Array{LayerBC}
    rCodeParams::Array{LayerBC}
    fitness::Real

    function Chromosome(m::Matrix, z::Array{LayerBC}, r::Array{LayerBC})
        if (z[1].min != 0 || z[1].max != 0)
            error("Range of depths for first layer must be (0,0)!")
        end
        if size(m)[1] != length(z) || size(m)[1] != length(r)
            error("Dimensions of zParams or rParams does not match number of layers")
        end

        new(m, size(m)[1], z, r, NaN)
    end
end

function createRandomModel(N::Integer, zParams::Array{LayerBC}, rParams::Array{LayerBC})
    if zParams[1].min != 0 || zParams[1].max != 0
        error("Range of depths for first layer must be (0,0)!")
    end

    if length(zParams) != N || length(rParams) != N
        error("Dimensions of zParams or rParams does not match number of layers")
    end

    model = zeros(N,2)
    model[1,2] = rand(rParams[1].min:10.0^-rParams[1].res:rParams[1].max)
    model[1,2] = signif(rand(rParams[1].min:rParams[1].max), rParams[1].res)
    for n = 2:N
        layerZParams = zParams[n]
        layerRParams = rParams[n]
        model[n,1] = signif(rand(layerZParams.min:layerZParams.max), layerZParams.res)
        model[n,2] = signif(rand(layerRParams.min:layerRParams.max), layerRParams.res)
    end

    Chromosome(model, zParams, rParams)
end

function calculateFitness!(c::Chromosome, data::Matrix)
    fs = data[:,1].^-1
    ρ, Φ = mt1d(c.model, fs)

    RMSρ = rms(data[:,2], ρ, data[:,4])
    RMSΦ = rms(data[:,3], Φ, data[:,5])

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
end

function crossover(a::Chromosome, b::Chromosome)
    const η = 2 # Distribution index
    cA = deepcopy(a)
    cB = deepcopy(b)
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

        # Apply significant figures
        cA.model[n,1] = signif(cA.model[n,1], layerZParams.res)
        cA.model[n,2] = signif(cA.model[n,2], layerRParams.res)
        cB.model[n,1] = signif(cB.model[n,1], layerZParams.res)
        cB.model[n,2] = signif(cB.model[n,2], layerRParams.res)
    end

    # Force first layer to have zero depth
    cA.model[1,1] = 0
    cB.model[1,1] = 0

    cA.model = sortrows(cA.model)
    cB.model = sortrows(cB.model)

    return (cA, cB)
end

# Computes β(u,a) and β(u,b) returning them in a tuple
# ## Arguments:
# - `p1` -- First parent value to crossover
# - `p2` -- Second parent value to crossover
# - `η`  -- Distribution parameter
# - `a`  -- Lower bound of search space
# - `b`  -- Upper bound of search space
#
# ## Returns:
# - A tuple containing β(u,a) and β(u,b) (in that order)
#
# ## Notes:
# Throws an error if a > b or if p1 > p2.
# Implementation is based on report by Deb, K. and Jain, H. (#2011017)
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

function mutate!(c::Chromosome, Pm::Real)
    for n in 1:c.N
        # Mutate depth
        if rand() < Pm
            # Do nothing to the depth if this is the top layer
            if n != 1
                zParams = c.zCodeParams[n]
                c.model[n,1] = signif(rand(zParams.min:zParams.max), zParams.res)
                sortrows(c.model)
            end
        end

        # Mutate resistivity
        if rand() < Pm
            rParams = c.rCodeParams[n]
            c.model[n,2] = signif(rand(rParams.min:rParams.max), rParams.res)
            sortrows(c.model)
        end
    end
end
