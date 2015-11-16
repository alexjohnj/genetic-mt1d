include("./mt1d.jl");
using MT1D
using PyPlot

function plotMT1DCompare(modelFile::AbstractString, dataFile::AbstractString)
    modelFileHandle = open(modelFile);
    model = readdlm(modelFileHandle);
    close(modelFileHandle);

    dataFileHandle = open(dataFile);
    data = readdlm(dataFileHandle);
    close(dataFileHandle);

    Ts     = data[:,1];
    fs     = Ts.^-1;
    ρData  = data[:,2];
    ΦData  = data[:,3];

    ρa, Φa = mt1d(model, fs);
    σρData = ρData * 0.12;
    σΦData = ΦData * 0.09;

    rmsΦ = rms(ΦData, Φa, σΦData);
    rmsρ = rms(ρData, ρa, σρData);
    rmsTotal = sqrt(rmsΦ^2 + rmsρ^2);
    titleStr = @sprintf("Total RMS: %.2f", rmsTotal);

    subplot(2,1,1);
    loglog(Ts, ρa);
    errorbar(Ts, ρData, yerr=σρData, fmt="x");
    xlabel("Period (S)");
    ylabel("Apparent Resistivity (Ohm-Metres)");
    title(titleStr);

    subplot(2,1,2);
    semilogx(Ts, Φa);
    errorbar(Ts, ΦData, yerr=σΦData, fmt="x");
    xlabel("Period (S)");
    ylabel("Phase (°)");

    legend(["Model", "Data"]);
end
