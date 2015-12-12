include("./MT1D.jl")

module MT1DGA
using MT1D

export LayerBC, Population, evolve!
"""
`LayerBC` defines a set of boundary conditions for a layer. One instance
represents either the resistivity or depth boundaries.

# Fields
- `min::Integer`: Lower boundary for the layer.
- `max::Integer`: Upper boundary for the layer.
"""
type LayerBC
    min::Integer
    max::Integer
end

include("./Model.jl")
include("./Population.jl")

end
