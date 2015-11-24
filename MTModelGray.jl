module MTModelGrayMod

export MTModelGray, encodeModel!, decodeModel!

type MTModelGray
    depthCode::AbstractString
    resCode::AbstractString
    N::Integer  # Number of Layers
    ρMin::Integer
    ρMax::Integer
    zMax::Integer
    Δz::Integer
    Δρ::Integer
    zNBits::Integer
    ρNBits::Integer

    function MTModelGray(dC,rC,N,ρMin,ρMax,zM,Δz,Δρ)
        zNBits = ceil(Integer, log2(zM / Δz))
        ρNBits = ceil(Integer, log2((ρMax - ρMin) / Δρ))
        new(dC,rC,N,ρMin,ρMax,zM,Δz,Δρ,zNBits,ρNBits)
    end
end

function encodeModel!(eM, model::Matrix)
    depthCode = ""
    resCode = ""
    for i = 1:eM.N
        zRelative = int2gray(Int(model[i,1] / eM.Δz))
        ρRelative = int2gray(Int((model[i,2] - eM.ρMin) / eM.Δρ))
        depthCode *=  bin(zRelative, eM.zNBits)
        resCode *= bin(ρRelative, eM.ρNBits)
    end
    eM.depthCode = depthCode
    eM.resCode = resCode
end

function decodeModel!(eM)
    model = zeros(Integer, (eM.N, 2));
    for i = 1:eM.N
        zLayerCode = "0b" * eM.depthCode[(i-1)*eM.zNBits + 1:i*eM.zNBits]
        resLayerCode = "0b" * eM.resCode[(i-1)*eM.ρNBits + 1:i * eM.ρNBits]
        zLayerVal = gray2int(parse(Int, zLayerCode)) * eM.Δz
        resLayerVal = eM.ρMin + (gray2int(parse(Int, resLayerCode)) * eM.Δρ)
        model[i,:] = [zLayerVal resLayerVal]
    end

    return model
end

function int2gray(n::Integer)
    n $ (n >> 1)
end

function gray2int(n::Integer)
    mask = n >> 1

    while mask != 0
        n = n $ mask
        mask = mask >> 1
    end
    return n
end
end
