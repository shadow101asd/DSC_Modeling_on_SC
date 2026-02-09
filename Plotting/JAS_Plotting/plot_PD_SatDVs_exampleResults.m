% Plotting the PD SatDV example plots

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

aF = 1*AU;
eF = 0.3;

maxN = 30;
Nws = 1:1:maxN;
NMs = 1:1:maxN;
frac = 1;
Npl = 5;
NL = 20;

[X,Y] = meshgrid(Nws,NMs);

%% Run

DVs = NaN(size(X));
tic
parfor i = 1:maxN
    Nw = Nws(i);
    for j = 1:maxN
        NM = NMs(j);
        SatDVs = compute_SIdeployment_SatDVs_2DF(Nw,NM,aF,eF,NL,Npl,frac);
        DVs(i,j) = max(SatDVs);
        % disp("Analyses completed: " + int2str((i-1)*A + j) + "/" + int2str(A^2))
    end
end
toc

%% Save Data

save("Data/plot_PD_SatDVs_exampleResults.mat", "DVs")

%% Plotting

w = 500;
h = 400;

n = 348729562;
f1 = figure(n);
h1 = imagesc(NMs,Nws,DVs);
set(h1, 'AlphaData', 1-isnan(DVs));
set(gca, ColorScale='linear', YDir='normal')
c1 = colorbar;
ylabel(c1,'max(\DeltaVs_{sats}) [km/s]', FontName="Times New Roman", FontSize=13);
colormap parula
xlabel('||N_M||')
ylabel('N_\omega')
set(gca, FontName="Times New Roman", FontSize=13)
set(f1, Position=[200 200 w h])


%% Saving 

filename = "Figures/PD_SatDVs_NwND_grid.eps";
exportgraphics(f1, filename);