"""
Description
===========

Sort a population in decreasing order of fitness.

Arguments
=========

- `P::Array{Model}`: The population to sort.

Side Effects
============

- Mutates the order of `P`.

Returns
=======

- `nothing`.
"""
function sortpop!(P::Array{Model})
    sort!(P, lt=(a,b) -> a.fitness < b.fitness)
    return
end

"""
Description
===========

Select a member of a population through a tournament selection.

Arguments
=========

- `P::Array{Model}`: Population to select a winner from.
- `K::Integer`: Number of competitors in tournament.
- `Ps::Real`: Probability of selection. 0 <= Ps <= 1.

Returns
=======

- `w::Model`: The winning model

Notes
=====

- **Throws an error** if `!(1 <= K <= length(P))`.
"""
function select_with_tournament(P::Array{Model}, K::Integer, Ps::Real)
    if !(1 <= K <= length(P))
        error("Tournament size can not be less than 1 or greater than the population size")
    end

    competitors = shuffle(P)[1:K]
    sortpop!(competitors)

    result = rand()
    for k in 0:K-2
        if result <= Ps*(1-Ps)^k
            return competitors[1]
        else
            competitors = competitors[2:end]
        end
    end

    return competitors[1]
end
