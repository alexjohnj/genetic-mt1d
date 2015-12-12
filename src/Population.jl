"""
Sort a population in decreasing order of fitness.

# Arguments

- `P::Array{Model}`: The population to sort.

# Side Effects

- Mutates the order of `P`.

# Returns

- `nothing`.
"""
function sortPopulation!(P::Array{Model})
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
"""
function performTournament(P::Array{Model}, K::Integer, Ps::Real)
    competitors = rand(P, K)
    sortPopulation!(competitors)

    tournamentResult = rand()
    for k in 0:K-2
        if tournamentResult <= Ps*(1-Ps)^k
            return competitors[1]
        else
            competitors = competitors[2:end]
        end
    end

    return competitors[1]
end
