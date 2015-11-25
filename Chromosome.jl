include("./mt1d.jl")
using MT1D
type LayerCodeParameters
    min::Integer
    max::Integer
    res::Integer
end

type Chromosome
    model::Matrix
    N::Integer
    zCodeParams::Array{LayerCodeParameters}
    rCodeParams::Array{LayerCodeParameters}
    fitness::Real

    function Chromosome(m::Matrix, z::Array{LayerCodeParameters}, r::Array{LayerCodeParameters})
        if (z[1].min != 0 || z[1].max != 0)
            error("Range of depths for first layer must be (0,0)!")
        end
        if size(m)[1] != length(z) || size(m)[1] != length(r)
            error("Dimensions of zParams or rParams does not match number of layers")
        end

        new(m, size(m)[1], z, r, NaN)
    end
end

function createRandomModel(N::Integer, zParams::Array{LayerCodeParameters}, rParams::Array{LayerCodeParameters})
    if zParams[1].min != 0 || zParams[1].max != 0
        error("Range of depths for first layer must be (0,0)!")
    end

    if length(zParams) != N || length(rParams) != N
        error("Dimensions of zParams or rParams does not match number of layers")
    end

    model = zeros(N,2)
    model[1,2] = rand(rParams[1].min:rParams[1].res:rParams[1].max)
    for n = 2:N
        layerZParams = zParams[n]
        layerRParams = rParams[n]
        model[n,1] = rand(layerZParams.min:layerZParams.res:layerZParams.max)
        model[n,2] = rand(layerRParams.min:layerRParams.res:layerRParams.max)
    end

    Chromosome(model, zParams, rParams)
end

function calculateFitness!(c::Chromosome)
    global data

    fs = data[:,1].^-1
    ρ, Φ = mt1d(c.model, fs)

    RMSρ = rms(data[:,2], ρ, data[:,4])^2
    RMSΦ = rms(data[:,3], Φ, data[:,5])^2

    c.fitness = sqrt(RMSρ^2 + RMSϕ^2)
end

function crossover(a::Chromosome, b::Chromosome)
    weight = rand()

    cA = deepcopy(a)
    cB = deepcopy(b)
    cA.model = weight * a.model + (1-weight) * b.model;
    cB.model = (1-weight) * a.model + (1-weight) * b.model;

    return(cA, cB)
end

function mutate!(c::Chromosome, Pm::Real)
    for n in 1:c.N
        # Mutate depth
        if rand() < Pm
            # Do nothing to the depth if this is the top layer
            if n != 1
                zParams = c.zCodeParams[n]
                c.model[n,1] = rand(zParams.min:zParams.res:zParams.max)
            end
        end

        # Mutate resistivity
        if rand() < Pm
            rParams = c.rCodeParams[n]
            c.model[n,2] = rand(rParams.min:rParams.res:rParams.max)
        end
    end
end
