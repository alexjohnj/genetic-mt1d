type Inversion
    data::Matrix
    popSize::Integer
    zBounds::Array{LayerBC}
    rBounds::Array{LayerBC}

    nE::Integer
    pM::Real
    pS::Real

    pop::Array{Model}
    gen::Integer

    function Inversion(data::Matrix, popSize::Integer, zBounds::Array{LayerBC},
                       rBounds::Array{LayerBC})
        Inversion(data, popSize, zBounds, rBounds, ceil(Integer, 0.05*popSize), 0.01, 0.65)
    end

    function Inversion(data::Matrix, popSize::Integer, zBounds::Array{LayerBC},
                       rBounds::Array{LayerBC}, nE::Integer, pM::Real, pS::Real)
        if length(zBounds) != length(rBounds)
            error("Length of zBounds != length of rBounds")
        end

        nLayers = length(zBounds)
        pop = vec(Array{Model}(popSize,1))
        for i in 1:popSize
            pop[i] = createRandomModel(nLayers, zBounds, rBounds)
            calculateFitness!(pop[i], data)
        end
        sortPopulation!(pop)
        new(data, popSize, zBounds, rBounds, nE, pM, pS, pop, 1)
    end
end

"""
"""
function evolve!(I::Inversion, ngen=100)
    @printf("Gen: %d\tBest RMS:%.2f\tSize: %d\n", I.gen, I.pop[1].fitness, length(I.pop))
    for i in 1:ngen
        createNextGeneration!(I)

        for idx in 1:I.popSize
            mutate!(I.pop[idx], I.pM)
            # TODO: Only recalculate fitness if mutation occurs
            calculateFitness!(I.pop[idx], I.data)
        end

        sortPopulation!(I.pop)
        I.gen += 1

        if 0.5 <= I.pop[1].fitness <= 1.5
            break
        end

        if I.gen % 100 == 0
            @printf("Gen: %d\tBest RMS:%.2f\tSize: %d\n", I.gen, I.pop[1].fitness, length(I.pop))
        end
    end
    return
end

"""
"""
function createNextGeneration!(I::Inversion)
    nextGen = I.pop[1:I.nE]

    # Create new children to fill the next generation
    while length(nextGen) < I.popSize
        parentA = performTournament(I.pop, 2, I.pS)
        parentB = performTournament(I.pop, 2, I.pS)

        childA, childB = crossover(parentA, parentB)
        push!(nextGen, childA, childB)
    end

    # At most, the next generation can have 1 child too many. Kill
    # this child if this is the case.
    if length(nextGen) > I.popSize
        I.pop = nextGen[1:end-1]
    else
        I.pop = nextGen
    end

    return
end
