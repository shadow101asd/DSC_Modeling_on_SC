% Evaluate Deliverable Dry Mass to a given 2D Flower Constellations,
% Using Different Insertion Strategies.

% Constants

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
Shuttle_specs = load("../Inputs/SE/SS_BIII.mat").Shuttle_specs;
Sat_Isp = 300; % s

% Parameters

NSats = 10;
M =  9;
W = -7;

NA = 2e3;
NE = 2e3;
aMIN = 0.1; % AU
aMAX = 6; % AU
eMIN = 1/(NE-1);
eMAX = 1-eMIN;
aFLOs = AU*linspace(aMIN, aMAX, NA);
eFLOs = linspace(eMIN, eMAX, NE);

%% Hohmann/2BI
% Setup data fields

Sat_masses = zeros([NA NE]);
TODs = zeros([NA NE]);
exitflags = zeros([NA NE]);


parfor i = 1:NA
    for j = 1:NE
        aFLO = aFLOs(i);
        eFLO = eFLOs(j);
        [m_dry_persat, TOD, exitflag] = getSatDryMasses2BI_F2D(aFLO,eFLO,M,W,muSu,NSats,Shuttle_specs.Isp,Sat_Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload)

        Sat_masses(i,j) = m_dry_persat;
        TODs(i,j) = TOD;
        exitflags(i,j) = exitflag;

        % Progress counter
        progress = int2str((i-1)*NE+j) + "/" + int2str(3*NA*NE) + " analyses run"
    end
end


%% Plotting Setup

[X,Y] = meshgrid(eFLOs, aFLOs/AU);
min_val = 1e2;
max_val = Shuttle_specs.maxPayload/NSats;
x = [min(X) max(X)];
y = [min(Y) max(Y)];


%% Hohmann Heatmaps

figure(121)
h1 = imagesc(x,y,Sat_masses);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Satellite Masses with NSats = " + num2str(NSats) + ", M = " + num2str(M) + ", W = " + num2str(W) + ". Hohmann Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("eFLO")
ylabel("aFLO [AU]")
ylim([0 ceil(max(aFLOs/AU))])

c1 = colorbar;
ylabel(c1,'Dry mass per satellite [kg]', FontSize=12);
colormap parula
fontsize(gcf, 14, 'points'); 

figure(122)
h1 = imagesc(x,y,TODs);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Total Deployment Time with NSats = " + num2str(NSats) + ", M = " + num2str(M) + ", W = " + num2str(W) + ". Hohmann Insertion")
set(gca,'YDir','normal');
xlabel("eFLO")
ylabel("aFLO [AU]")
ylim([0 ceil(max(aFLOs/AU))])

c2 = colorbar;
ylabel(c2,'Total Deployment Time [Earth years]', FontSize=12);
colormap parula
fontsize(gcf, 14, 'points'); 

figure(123)
plot_F2D_Static(123, 1.25*AU, 0.35, M, W, muSu, "plot", false)