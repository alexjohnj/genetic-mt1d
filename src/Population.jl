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
Select a member of the population using a tournament selection. A
random number is generated. If it is less than the probability of
selection, the fittest competitor wins. Otherwise, the fittest is
discarded and the tournament is carried out between the remaining
weaker candidates. If no body wins the tournament, the weakest
candidate wins by default.

## Key Word Arguments
- P  -- Population to select competitors from
- K  -- Tournament size
- Ps -- Probability of selection

## Returns:
A winning `Model`
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
