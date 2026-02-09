% Evaluate Deliverable Dry Mass to a given 2D Flower Constellations,
% Using the Phasing Insertion Strategy

% Add paths
addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts
addpath("SatM_funcs/") % Helpers

% Constants

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
Shuttle_specs = load("../../Inputs/SE/SS_BIII.mat").Shuttle_specs;
Sat_Isp = 300; % s

% Parameters

Npl = 10;
NL = 10;
f = 1;
Nw =  3;
NM =  2;

NA = 100;
NE = 100;
aMIN = 0.1; % AU
aMAX = 5.5; % AU
eMIN = 1/(NE-1);
eMAX = 1-eMIN;
aFs = AU*linspace(aMIN, aMAX, NA);
eFs = linspace(eMIN, eMAX, NE);

%% Preprocessing for fixed aF = AU

resolution = 500;

Sat_DVs = NaN(NE, Npl);
DPs_opt = NaN(NE, Npl);
Sat_TOFs = NaN(NE, Npl);
tic
parfor i = 1:NE
    [Sat_DVs_i, DPs_opt_i] = compute_SIdeployment_SatDVs_2DF(Nw,NM,AU,eFs(i),NL,Npl,f,resolution);
    Sat_DVs(i,:) = Sat_DVs_i*1e3; % m/s
    DPs_opt(i,:) = DPs_opt_i;

    Sat_TOFs_i = NaN(1,Npl);
    for s = 1:Npl
        [~, ~, ~, ~, details] = computeOptimalPhasingDV_2DF_extended(Nw, -NM, eFs(i), DPs_opt_i(s));
        Sat_TOFs_i(s) = details.t2;
    end
    Sat_TOFs(i,:) = Sat_TOFs_i;
end
toc

%% Saving preprocessing

filename_PP = "./Data/PhasingD_PP_3.mat";
save(filename_PP, "Sat_DVs","DPs_opt", "Sat_TOFs")

%% Phasing Insertion
% Setup data fields

Sat_masses = zeros([NA NE]);
TODs = zeros([NA NE]);
exitflags = zeros([NA NE]);

load(filename_PP, "Sat_DVs","DPs_opt", "Sat_TOFs")

tic
parfor i = 1:NA
    aFLO = aFs(i);
    
    for j = 1:NE
        eFLO = eFs(j);

        SDVs = Sat_DVs(j,:) * (aFLO/AU)^(-1/2); % scale appropriately
        SDPs = DPs_opt(j,:);
        STOFs = Sat_TOFs(j,:) * (aFLO/AU)^(3/2); % scale appropriately

        [m_dry_persat, TOD, exitflag] = getSatDryMassesPhasingD_F2D(aFLO,eFLO,NM,Nw,muSu,NL,Npl,f, ...
            Shuttle_specs.Isp,Sat_Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, ...
            SDVs, SDPs, STOFs)

        Sat_masses(i,j) = m_dry_persat;
        TODs(i,j) = TOD;
        exitflags(i,j) = exitflag;

        % Progress counter
        progress = int2str((i-1)*NE+j) + "/" + int2str(NA*NE) + " analyses run"
    end
end
toc

%% Plotting

[X,Y] = meshgrid(eFs, aFs/AU);
min_val = 1e2;
max_val = Shuttle_specs.maxPayload/Npl;
x = [min(X) max(X)];
y = [min(Y) max(Y)];

w = 500;
h = 300;

% Hohmann Heatmaps

n = 177725732;
f1 = figure(n);
h1 = imagesc(x,y,Sat_masses);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
% set(gca,'YDir','normal');
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("e_F")
ylabel("a_F [AU]")
ylim([0 ceil(max(aFs/AU))])

c1 = colorbar;
ylabel(c1,'Dry mass per satellite [kg]', FontSize=12, FontName="Times New Roman");
colormap parula
fontsize(f1, 14, 'points');
set(f1, Position=[0 h w h]);
set(gca, FontName="Times New Roman")

f2 = figure(n+1);
h1 = imagesc(x,y,TODs);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
set(gca,'YDir','normal');
xlabel("e_F")
ylabel("a_F [AU]")
ylim([0 ceil(max(aFs/AU))])

c2 = colorbar;
ylabel(c2,'Total Deployment Time [years]', FontSize=12, FontName="Times New Roman");
colormap parula
fontsize(f2, 14, 'points');
set(f2, Position=[w h w h]); 
set(gca, FontName="Times New Roman")

%% Saving

filename1 = "./Figures/SatMass_HeatMap_PhasingDeployment.eps";
exportgraphics(f1, filename1);

% filename2 = "./Figures/TOD_HeatMap_ApoDeployment.eps";
% exportgraphics(filename2, f2);

%% Saving data

filename_data = "./Data/SatPhasingD.mat";
save(filename_data, "Sat_masses", "TODs", "exitflags")