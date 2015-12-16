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
- `popsize::Integer`: The target size of the population.
- `depth_lims::Array{LayerBC}`: Array of boundaries for each layer's depth.
- `res_lims::Array{LayerBC}`: Array of boundaries for each layer's resistivity.
- `nelitists::Integer`: Number of elitist clones to make per generation.
- `pmutation::Real`: Probability of mutation. 0 <= Pm <= 1.
- `pselection::Real`: Probability of selection. 0 <= Ps <= 1.
- `tourn_size::Integer`: Number of competitors in tournament selection
- `η::Integer`: Distribution index for SBX. Must be non-negative.
- `pop::Array{Model}`: Population for the genetic algorithm. Sorted by fitness.
- `gen::Integer`: Current generation. Starts at 1.

Usage
=====

- `I = Inversion(data, popsize, depth_lims, res_lims)`
- `I = Inversion(data, popsize, depth_lims, res_lims, nelitists, pmutation, pselection, tourn_size, η)`


Notes
=====

- If `nelitists` is unspecified, it defaults to 5% of the population.
- If `pmutation` is unspecified, it defaults to 0.01.
- If `pselection` is unspecified, it defaults to 0.65.
- If `tourn_size` is unspecified, it defaults to 2.
- If `η` is unspecified, it defaults to 2.
"""
type Inversion
    data::Matrix
    popsize::Integer
    depth_lims::Array{LayerBC}
    res_lims::Array{LayerBC}

    nelitists::Integer
    pmutation::Real
    pselection::Real
    tourn_size::Integer

    η::Integer

    pop::Array{Model}
    gen::Int

    function Inversion(data::Matrix, popsize::Integer, depth_lims::Array{LayerBC},
                       res_lims::Array{LayerBC})
        Inversion(data, popsize, depth_lims, res_lims, ceil(Integer, 0.05 * popsize), 0.01, 0.65, 2, 2)
    end

    function Inversion(data::Matrix, popsize::Integer, depth_lims::Array{LayerBC},
                       res_lims::Array{LayerBC}, nelitists::Integer, pmutation::Real, pselection::Real,
                       tourn_size::Int, η::Integer)
        if length(depth_lims) != length(res_lims)
            error("Length of depth_lims != length of res_lims")
        end

        nLayers = length(depth_lims)
        pop = vec(Array{Model}(popsize,1))
        for i in 1:popsize
            pop[i] = rand_model(nLayers, depth_lims, res_lims)
            calcfitness!(pop[i], data)
        end
        sortpop!(pop)
        new(data, popsize, depth_lims, res_lims, nelitists, pmutation, pselection, tourn_size, η, pop, 1)
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
- `verbose::Bool=false`: Print the best RMS misfit every 100 generations.

Side Effects
============

- Mutates the `pop` and `gen` fields of `I`.

Returns
=======

- `nothing`
"""
function evolve!(I::Inversion, ngen=100; verbose=false)
    if verbose
        @printf("Gen: %d\tBest RMS:%.2f\tSize: %d\n", I.gen, I.pop[1].fitness, length(I.pop))
    end
    for i in 1:ngen
        make_new_generation!(I)

        for idx in 1:I.popsize
            mutate!(I.pop[idx], I.pmutation)
            calcfitness!(I.pop[idx], I.data)
        end

        sortpop!(I.pop)
        I.gen += 1

        if 0.5 <= I.pop[1].fitness <= 1.5
            break
        end

        if verbose && I.gen % 100 == 0
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
function make_new_generation!(I::Inversion)
    nextgen = I.pop[1:I.nelitists]

    # Create new children to fill the next generation
    while length(nextgen) < I.popsize
        p1 = select_with_tournament(I.pop, I.tourn_size, I.pselection)
        p2 = select_with_tournament(I.pop, I.tourn_size, I.pselection)

        childA, childB = crossover(p1, p2, I.η)
        calcfitness!(childA, I.data)
        calcfitness!(childB, I.data)
        push!(nextgen, childA, childB)
    end

    # At most, the next generation can have 1 child too many. Kill
    # this child if this is the case.
    if length(nextgen) > I.popsize
        I.pop = nextgen[1:end-1]
    else
        I.pop = nextgen
    end

    return
end
