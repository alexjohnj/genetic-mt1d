include("./MT1D.jl")

module MT1DGA
using MT1D

export LayerBC, Population, evolve!
type LayerBC
    min::Integer
    max::Integer
end

include("./Chromosome.jl")
include("./Population.jl")

end
