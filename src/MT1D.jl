module MT1D
export mt1d
export rms

"""
MT1D calculates the apparent resistivity and phase for an `N` layer 1D
resistivity model.

## Arguments:
- `model` -- An `N`×2 matrix with depth to the top of each layer in the 1st
    column and resistivity in the second.
- `fs` -- A vector of frequencies to evaluate the model at.

## Returns:
- `ρa` -- A vector of apparent surface resistivities computed from the model.
- `Φ`  -- A vector of phases computed from the model.
"""
function mt1d(model::Matrix, fs::Vector)
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
rms finds the root mean squared misfit between a set of observations and a
model while accounting for uncertainty.

## Arguments:
- `obs`   -- A vector f observations.
- `model` -- A vector of models.
- `ϵ`     -- A vector of uncertainties.

## Returns:
- A real number indicating the RMS value.
"""
function rms(obs::Vector, model::Vector, ϵ::Vector)
    sqrt(sum(((obs - model) ./ ϵ).^2) / (2 * length(obs)))
end
end
