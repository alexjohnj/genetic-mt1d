"""
Description
===========

An instance of `Inversion` holds the current state of an inversion. This
includes the population, the current generation, the boundaries of the inversion
and the inversion parameters.

Fields
======

- `data::Matrix`: The data being inverted. Columns should consist of:
    1. Period
    2. Apparent resistivity
    3. Phase
    4. Apparent resistivity uncertainty
    5. Phase uncertainty
- `popSize::Integer`: The target size of the population.
- `zBounds::Array{LayerBC}`: Array of boundaries for each layer's depth.
- `rBounds::Array{LayerBC}`: Array of boundaries for each layer's resistivity.
- `nE::Integer`: Number of elitist clones to make per generation.
- `pM::Real`: Probability of mutation. 0 <= Pm <= 1.
- `pS::Real`: Probability of selection. 0 <= Ps <= 1.
- `tSize::Integer`: Number of competitors in tournament selection
- `pop::Array{Model}`: Population for the genetic algorithm. Sorted by fitness.
- `gen::Integer`: Current generation. Starts at 1.

Usage
=====

- `I = Inversion(data, popSize, zBounds, rBounds)`
- `I = Inversion(data, popSize, zBounds, rBounds, nE, pM, pS, tSize)`


Notes
=====

- If `nE` is unspecified, it defaults to 5% of the population.
- If `pM` is unspecified, it defaults to 0.01.
- If `pS` is unspecified, it defaults to 0.65.
- If `tSize` is unspecified, it defaults to 2.
"""
type Inversion
    data::Matrix
    popSize::Integer
    zBounds::Array{LayerBC}
    rBounds::Array{LayerBC}

    nE::Integer
    pM::Real
    pS::Real
    tSize::Integer

    pop::Array{Model}
    gen::Int

    function Inversion(data::Matrix, popSize::Integer, zBounds::Array{LayerBC},
                       rBounds::Array{LayerBC})
        Inversion(data, popSize, zBounds, rBounds, ceil(Integer, 0.05*popSize), 0.01, 0.65, 2)
    end

    function Inversion(data::Matrix, popSize::Integer, zBounds::Array{LayerBC},
                       rBounds::Array{LayerBC}, nE::Integer, pM::Real, pS::Real,
                       tSize::Int)
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
        new(data, popSize, zBounds, rBounds, nE, pM, pS, tSize, pop, 1)
    end
end

"""
Description
===========

Evolve an instance of `Inversion` through a specified number of
generations. This involves carrying out all the steps of the genetic
algorithm. Evolution will be terminated if the fitness falls within the stopping
conditions.

Arguments
=========

- `I::Inversion`: Inversion instance to evolve.
- `ngen::Integer=100`: Number of generations to advance `I` by.

Side Effects
============

- Mutates the `pop` and `gen` fields of `I`.

Returns
=======

- `nothing`
"""
function evolve!(I::Inversion, ngen=100)
    @printf("Gen: %d\tBest RMS:%.2f\tSize: %d\n", I.gen, I.pop[1].fitness, length(I.pop))
    for i in 1:ngen
        createNextGeneration!(I)

        for idx in 1:I.popSize
            mutate!(I.pop[idx], I.pM)
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
Description
===========

Create the next generation for an inversion by making elitist clones, selecting
parents and mating the selections.

Arguments
=========

- `I::Inversion`: Inversion instance to create the next generation for.

Side Effects
============

- Mutates the `pop` field of `I`.

Returns
=======

- `nothing`
"""
function createNextGeneration!(I::Inversion)
    nextGen = I.pop[1:I.nE]

    # Create new children to fill the next generation
    while length(nextGen) < I.popSize
        parentA = performTournament(I.pop, I.tSize, I.pS)
        parentB = performTournament(I.pop, I.tSize, I.pS)

        childA, childB = crossover(parentA, parentB)
        calculateFitness!(childA, I.data)
        calculateFitness!(childB, I.data)
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
