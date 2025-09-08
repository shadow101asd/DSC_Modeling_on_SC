% Evaluate Deliverable Dry Mass to a range of 2D Flower Constellations,
% Using Different Insertion Strategies. Adapt aFLO, eFLO, to be synodic
% with a planet of specified semi-major axis

% Constants

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
Shuttle_specs = load("../Inputs/SE/SS_BIII.mat").Shuttle_specs;
Sat_Isp = 300; % s

% Parameters

NSats = 10;
aPL = 1.0*AU; % Earth: 1.000 AU. Mars: 1.524 AU. Venus: 0.723 AU.

maxF = 30;
maxW = 30;
Fs = 1:maxF;
Ws = 1:maxW;

%% Hohmann/2BI
% Setup data fields

Sat_masses = zeros([maxF maxW]);
TODs = zeros([maxF maxW]);
exitflags = zeros([maxF maxW]);
aFLOs = zeros([maxF maxW]);
eFLOs = zeros([maxF maxW]);


parfor i = Fs
    for j = Ws
        GCD = gcd(i,j);
        F = i/GCD;
        W = j/GCD;

        aFLO = aPL * abs(W/F)^(2/3);
        eFLO = abs(aPL/aFLO - 1);

        if eFLO >= 1 || eFLO < 0
            Sat_masses(i,j) = NaN;
            TODs(i,j) = NaN;
            exitflags(i,j) = 0;
            aFLOs(i,j) = NaN;
            eFLOs(i,j) = NaN;
        else
            [m_dry_persat, TOD, exitflag] = getSatDryMasses2BI_F2D(aFLO,eFLO,F,-W,muSu,NSats,Shuttle_specs.Isp,Sat_Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload)

            Sat_masses(i,j) = m_dry_persat;
            TODs(i,j) = TOD;
            exitflags(i,j) = exitflag;
            aFLOs(i,j) = aFLO;
            eFLOs(i,j) = eFLO;
        end

        

        % Progress counter
        progress = int2str((i-1)*maxW+j) + "/" + int2str(3*maxF*maxW) + " analyses run"
    end
end


%% Plotting Setup

[X,Y] = meshgrid(Ws, Fs);
min_val = 1e2;
max_val = Shuttle_specs.maxPayload/NSats;
x = [min(X) max(X)];
y = [min(Y) max(Y)];


%% Hohmann Heatmaps

figure(1121)
h1 = imagesc(x,y,Sat_masses);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Satellite Masses with NSats = " + num2str(NSats) + ", Hohmann Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("||W||")
ylabel("||F||")
xlim([1 maxW])
ylim([1 maxF])
fontsize(gcf, 14, 'points'); 

c1 = colorbar;
ylabel(c1,'Dry mass per satellite [kg]', FontSize=15);
colormap parula

figure(1122)
h1 = imagesc(x,y,TODs);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Total Deployment Time with NSats = " + num2str(NSats) + ", Hohmann Insertion")
set(gca,'YDir','normal');
xlabel("||W||")
ylabel("||F||")
xlim([1 maxW])
ylim([1 maxF])
fontsize(gcf, 14, 'points'); 

c2 = colorbar;
ylabel(c2,'Total Deployment Time [Earth years]', FontSize=15);
colormap parula

figure(1123)
h1 = imagesc(x,y,aFLOs/AU);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Selected aFLO [AU]")
set(gca,'YDir','normal');
xlabel("||W||")
ylabel("||F||")
xlim([1 maxW])
ylim([1 maxF])
fontsize(gcf, 14, 'points'); 

c3 = colorbar;
ylabel(c3,'Selected aFLO [AU]', FontSize=15);
colormap parula

figure(1124)
h1 = imagesc(x,y,eFLOs);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("2DFlower Selected eFLO")
set(gca,'YDir','normal');
xlabel("||W||")
ylabel("||F||")
xlim([1 maxW])
ylim([1 maxF])
fontsize(gcf, 14, 'points'); 

c4 = colorbar;
ylabel(c4,'Selected eFLO', FontSize=15);
colormap parula