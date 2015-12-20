module MT1D
export calc_response, calc_rms

"""
Description
===========

Forward model the apparent resistivity and phase response for a model.

Arguments
=========

- `model::Matrix`: Matrix representing the model. Each row is a layer. Columns are:
    1. Depth to top of layer (m)
    2. Resistivity of layer (Ωm)
- `fs::Vector`: Frequencies to forward model the response at.

Returns
=======

- `Tuple{ρ::Vector,Φ::Vector}`: Forward modelled apparent resistivity and phase.
"""
function calc_response{T<:AbstractFloat}(model::Matrix{T}, fs::Vector{T})
    const μ = 4E-7 * π
    nf = length(fs)
    nlayers = length(model[:,1])

    ρa = zeros(nf)
    Φ = zeros(nf)

    for fidx in 1:nf
        k = sqrt(-im * 2π * fs[fidx] * μ / model[end, 2])
        C = 1 / k
        for n in nlayers-1:-1:1
            k = sqrt(-im * 2π * fs[fidx] * μ / model[n, 2])
            Δz = model[n+1, 1] - model[n, 1]
            C = (C * k + tanh(k * Δz)) / (C * k^2 * tanh(k * Δz) + k)
        end
        Z = 2π * fs[fidx] * μ * C
        ρa[fidx] = abs(Z)^2 / (2π * fs[fidx] * μ)
        Φ[fidx] = 90 - angle(Z) * 180 / π
    end

    return (ρa, Φ);
end

"""
Description
===========

Calculate the root mean squared misfit between a set of observations and
modelled responses.

Arguments
=========

- `obs::Vector`: Observation points
- `model::Vector`: Modelled responses
- `ϵ::Vector`: Uncertainties for `obs`

Returns
=======

- `r::Real`: Calculated RMS misfit
"""
function calc_rms(obs::Vector, model::Vector, ϵ::Vector)
    sqrt(sum(((obs - model) ./ ϵ).^2) / (2 * length(obs)))
end
end
