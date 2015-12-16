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
function calc_response(model::Matrix, fs::Vector)
    const μ0 = 4E-7 * π;

    zs = model[:,1];
    ρs = model[:,2];
    N = length(model[:,1]);

    # Make sure zs is a row vector and ρs is a column vector
    fs = reshape(fs, (1, length(fs)));
    ρs = vec(ρs);

    # N x length(fs) matrix with k(row,col) corresponding to k[N,fsIdx]
    k = sqrt(-im * 8E-7 * π^2 * (ρs.^-1) * fs);
    C = 1 ./ k[end,:];
    Δz = zs[2:end] - zs[1:end-1];

    for jj = N-1:-1:1
        C = (C .* k[jj, :] + tanh(Δz[jj] * k[jj,:])) ./ (k[jj,:].^2 .* C .* tanh(Δz[jj] * k[jj,:]) + k[jj,:]);
    end

    # Calculate surface impedance, resistivity and phase
    Z  = μ0 * 2π * fs .* C[1,:];
    ρa = (2π * fs * μ0).^-1 .* abs(Z).^2;
    Φ  = 90 - angle(Z) * 180 / π;

    return (vec(ρa), vec(Φ));
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
