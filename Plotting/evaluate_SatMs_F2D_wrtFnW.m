% Evaluate Deliverable Dry Mass to a range of 2D Flower Constellations,
% Using Different Insertion Strategies. (for a given aFLO, eFLO)

% Constants

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
Shuttle_specs = load("../Inputs/SE/SS_BIII.mat").Shuttle_specs;
Sat_Isp = 300; % s

% Parameters

NSats = 10;
aFLO = 1.25*AU;
eFLO = 0.25;

maxM = 30;
maxW = 30;
Ms = 1:maxM;
Ws = 1:maxW;

%% Hohmann/2BI
% Setup data fields

Sat_masses = zeros([maxM maxW]);
TODs = zeros([maxM maxW]);
exitflags = zeros([maxM maxW]);


parfor i = Ms
    for j = Ws
        GCD = gcd(i,j);
        M = i/GCD;
        W = j/GCD;


        [m_dry_persat, TOD, exitflag] = getSatDryMasses2BI_F2D(aFLO,eFLO,M,-W,muSu,NSats,Shuttle_specs.Isp,Sat_Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload)

        Sat_masses(i,j) = m_dry_persat;
        TODs(i,j) = TOD;
        exitflags(i,j) = exitflag;

        % Progress counter
        progress = int2str((i-1)*maxW+j) + "/" + int2str(3*maxM*maxW) + " analyses run"
    end
end


%% Plotting Setup

[X,Y] = meshgrid(Ws, Ms);
min_val = 1e2;
max_val = Shuttle_specs.maxPayload/NSats;
x = [min(X) max(X)];
y = [min(Y) max(Y)];


%% Hohmann Heatmaps

figure(121)
h1 = imagesc(x,y,Sat_masses);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Satellite Masses with NSats = " + num2str(NSats) + ", Hohmann Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("||W||")
ylabel("||M||")
xlim([1 maxW])
ylim([1 maxM])
fontsize(gcf, 14, 'points'); 

c1 = colorbar;
ylabel(c1,'Dry mass per satellite [kg]', FontSize=15);
colormap parula

figure(122)
h1 = imagesc(x,y,TODs);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Total Deployment Time with NSats = " + num2str(NSats) + ", Hohmann Insertion")
set(gca,'YDir','normal');
xlabel("||W||")
ylabel("||M||")
xlim([1 maxW])
ylim([1 maxM])
fontsize(gcf, 14, 'points'); 

c2 = colorbar;
ylabel(c2,'Total Deployment Time [Earth years]', FontSize=15);
colormap parula