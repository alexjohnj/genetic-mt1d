type Population
    Q::Integer
    cs::Array{Chromosome}
    data::Matrix
    zBounds::Array{LayerBC}
    rBounds::Array{LayerBC}
    generation::Integer

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
        sort!(P, lt=(a,b) -> a.fitness < b.fitness)
        new(Q, P, data, zBounds, rBounds, 1)
    end
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
A winning `Chromosome`
"""
function performTournament(P::Population, K::Integer, Ps::Real)
    competitors = rand(P.cs, K)
    sort!(competitors, lt=(a,b) -> a.fitness < b.fitness)

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

"""
Creates the next generation by selecting a percentage of the
population to survive and filling the remaining population with
children. The children's parents are selected using a tournament.

## Arguments:
- P  -- The population to advance
- Pe -- Elitism rate

## Returns:

- Nothing
"""
function createNextGeneration!(P::Population, Pe::Real)
    if !(0 < Pe <= 1)
        error("Elitism parameter was not in the range 0<Pe<=1")
    end
    nToCopy = floor(Int, P.Q * Pe)
    nextGen = P.cs[1:nToCopy]

    # Create new children to fill the next generation
    while length(nextGen) < P.Q
        parentA = performTournament(P, 2, 0.6)
        parentB = performTournament(P, 2, 0.6)

        childA, childB = crossover(parentA, parentB)
        push!(nextGen, childA, childB)
    end

    # At most, the next generation can have 1 child too many. Kill
    # this child if this is the case.
    if length(nextGen) > P.Q
        P.cs = nextGen[1:end-1]
    else
        P.cs = nextGen
    end
end

function evolve!(P::Population, maxGen=1000)
    @printf("Generation: %d\tBest RMS: %.2f\tSize: %d\n", P.generation, P.cs[1].fitness, length(P.cs))
    while P.generation < maxGen
        createNextGeneration!(P, 0.1)
        for i = 1:P.Q
            mutate!(P.cs[i], 0.01)
            calculateFitness!(P.cs[i], P.data)
        end

        sort!(P.cs, lt=(a,b) -> a.fitness < b.fitness)
        P.generation += 1
        if 0.5 <= P.cs[1].fitness <= 1.5
            break
        end

        if P.generation % 100 == 0
            @printf("Generation: %d\tBest RMS: %.2f\tSize: %d\n", P.generation, P.cs[1].fitness, length(P.cs))
        end
    end
    println("*"^80)
    @printf("Generation: %d\tBest RMS: %.2f\tSize: %d\n", P.generation, P.cs[1].fitness, length(P.cs))
end
