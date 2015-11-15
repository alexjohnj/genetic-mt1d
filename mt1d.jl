const μ0 = 4E-7 * π;

function MT1D(model::Matrix{Float64}, fs::Vector{Float64})
    # Model = Nx2 matrix with depths in column 1 and resistivities in column 2
    # fs = Vector of frequencies

    zs = model[:,1];
    ρs = model[:,2];
    N = length(model[:,1]);

    # Make sure zs is a row vector and ρs is a column vector
    fs = reshape(fs, (1, length(fs)));
    ρs = vec(ρs);

    # N x len(fs) matrix with k(row,col) corresponding to k(N,fsIdx)
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
