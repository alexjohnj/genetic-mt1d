include("./Chromosome.jl")

type Population
    Q::Integer
    cs::Array{Chromosome}
    data::Matrix
    zBounds::Array{LayerBC}
    rBounds::Array{LayerBC}

    function Population(Q::Integer, data::Matrix, zBounds::Array{LayerBC}, rBounds::Array{LayerBC})
        if length(zBounds) != length(rBounds)
            error("Length of BC arrays did not match!")
        end
        # Number of layers
        N = size(zBounds)[1]
        P = []
        for i in 1:Q
            push!(P, createRandomModel(N, zBounds, rBounds))
            calculateFitness!(P[i], data)
        end
        new(Q, P, data, zBounds, rBounds)
    end
end
