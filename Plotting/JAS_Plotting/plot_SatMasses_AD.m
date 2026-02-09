% Evaluate Deliverable Dry Mass to a given 2D Flower Constellations,
% Using the Apoapsis Insertion Strategy

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

NA = 2e3;
NE = 2e3;
aMIN = 0.1; % AU
aMAX = 5.5; % AU
eMIN = 1/(NE-1);
eMAX = 1-eMIN;
aFs = AU*linspace(aMIN, aMAX, NA);
eFs = linspace(eMIN, eMAX, NE);

%% 2BI (Apoapsis)
% Setup data fields

Sat_masses = zeros([NA NE]);
TODs = zeros([NA NE]);
exitflags = zeros([NA NE]);


parfor i = 1:NA
    for j = 1:NE
        aFLO = aFs(i);
        eFLO = eFs(j);
        [m_dry_persat, TOD, exitflag] = getSatDryMassesAD_F2D(aFLO,eFLO,NM,Nw,muSu,NL,Npl,f,Shuttle_specs.Isp,Sat_Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload)

        Sat_masses(i,j) = m_dry_persat;
        TODs(i,j) = TOD;
        exitflags(i,j) = exitflag;

        % Progress counter
        progress = int2str((i-1)*NE+j) + "/" + int2str(3*NA*NE) + " analyses run"
    end
end


%% Plotting

[X,Y] = meshgrid(eFs, aFs/AU);
min_val = 1e2;
max_val = Shuttle_specs.maxPayload/Npl;
x = [min(X) max(X)];
y = [min(Y) max(Y)];

w = 500;
h = 300;

% Hohmann Heatmaps

n = 12352532;
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

filename1 = "./Figures/SatMass_HeatMap_ApoDeployment.eps";
exportgraphics(f1, filename1);

% filename2 = "./Figures/TOD_HeatMap_ApoDeployment.eps";
% exportgraphics(filename2, f2);