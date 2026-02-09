% Plotting the PD SatDV example plots

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

eps = 1e-2;
aF = 1*AU;
E = 100;
eFs = linspace(eps, 1-eps, E);

F = 100;
fracs = linspace(eps,1, F);

Nw = 3;
NM = 2;

Npl = 5;
NL = 1;


%% Run

DVs = NaN([F,E]);
tic
for i = 1:F
    disp("Currently on frac element "+ int2str(i)+"/"+int2str(F))
    f = fracs(i);
    parfor j = 1:E
        eF = eFs(j);
        SatDVs = compute_SIdeployment_SatDVs_2DF(Nw,NM,aF,eF,NL,Npl,f);
        DVs(i,j) = max(SatDVs);
        % disp("Analyses completed: " + int2str((i-1)*A + j) + "/" + int2str(A^2))
    end
end
toc

%% Save Data

save("Data/plot_PD_SatDVs_eFsvsfs_2.mat", "DVs")

%% Plotting

w = 500;
h = 400;

n = 348729562;
f1 = figure(n);
h1 = imagesc(eFs,fracs,DVs);
set(h1, 'AlphaData', 1-isnan(DVs));
set(gca, ColorScale='linear', YDir='normal')
c1 = colorbar;
ylabel(c1,'max(\DeltaVs_{sats}) [km/s]', FontName="Times New Roman", FontSize=13);
colormap parula
xlabel('e_F')
ylabel('f')
set(gca, FontName="Times New Roman", FontSize=13)
set(f1, Position=[200 200 w h])


%% Saving 

filename = "Figures/PD_SatDVs_eFvsf_grid_2.eps";
exportgraphics(f1, filename);