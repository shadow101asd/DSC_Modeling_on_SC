% Evaluate 2BI Dry Masses

% Assumed performance parameters

StarshipBIII_dryM =  100000; % kg
StarshipBIII_wetM = 2300*10^3+StarshipBIII_dryM; % kg
raptorVac_Isp = 382; % s
sat_Isp = 300; % s
maxSSPayload2LEO = 200*10^3; % kg

% Constants

muSu = 1.327124400419393e+11;
AU = 149600000; % 1 AU in km

% Analysis parameters

NSats = 10;
NA = 500;
NE = 500;
aMOGs = AU*linspace(0.4, 6, NA);
eMOGs = linspace(0.001, 0.99, NE);

%% Hohmann
% Setup data fields

Sat_masses = zeros([NA NE]);
Shuttle_FEs = zeros([NA NE]);
Sat_fuels = zeros([NA NE]);
leftoverShuttle_FMs = zeros([NA NE]);


for i = 1:NA
    for j = 1:NE
        aMOG = aMOGs(i);
        eMOG = eMOGs(j);
        [Sats_MOGM, Shuttle_FE, Sat_FE, leftoverShuttle_FM] = ...
            getSatDryMasses2BI(aMOG,eMOG,muSu,NSats,raptorVac_Isp,sat_Isp,StarshipBIII_wetM,StarshipBIII_dryM,maxSSPayload2LEO);
        Sat_masses(i,j) = Sats_MOGM(1);
        % Shuttle_FEs(i,j) = Shuttle_FE;
        % Sat_fuels(i,j) = Sat_FE;
        % leftoverShuttle_FMs(i,j) = leftoverShuttle_FM;

        % Progress counter
        progress = int2str((i-1)*NE+j) + "/" + int2str(3*NA*NE) + " analyses run"
    end
end

%% 1BI
% Setup data fields

Sat_masses_1BI = zeros([NA NE]);
Shuttle_FEs_1BI = zeros([NA NE]);
Sat_fuels_1BI = zeros([NA NE]);
leftoverShuttle_FMs_1BI = zeros([NA NE]);

for i = 1:NA
    for j = 1:NE
        aMOG = aMOGs(i);
        eMOG = eMOGs(j);
        [Sats_MOGM, Shuttle_FE, Sat_FE, leftoverShuttle_FM] = ...
            getSatDryMasses1BI(aMOG,eMOG,muSu,NSats,raptorVac_Isp,sat_Isp,StarshipBIII_wetM,StarshipBIII_dryM,maxSSPayload2LEO);
        Sat_masses_1BI(i,j) = Sats_MOGM(1);
        % Shuttle_FEs_1BI(i,j) = Shuttle_FE;
        % Sat_fuels_1BI(i,j) = Sat_FE;
        % leftoverShuttle_FMs_1BI(i,j) = leftoverShuttle_FM;

        % Progress counter
        progress = int2str((i-1)*NE+j + NA*NE) + "/" + int2str(3*NA*NE) + " analyses run"
    end
end

%% Spreadout Insertion ("SI")
% Setup data fields

Sat_masses_SI = zeros([NA NE]);
Shuttle_FEs_SI = zeros([NA NE]);
Sat_fuels_SI = zeros([NA NE]);
leftoverShuttle_FMs_SI = zeros([NA NE]);

for i = 1:NA
    for j = 1:NE
        aMOG = aMOGs(i);
        eMOG = eMOGs(j);
        [Sats_MOGM, Shuttle_FE, Sat_FE, leftoverShuttle_FM] = ...
            getSatDryMassesSI(aMOG,eMOG,muSu,NSats,raptorVac_Isp,sat_Isp,StarshipBIII_wetM,StarshipBIII_dryM,maxSSPayload2LEO);
        Sat_masses_SI(i,j) = Sats_MOGM(1);
        % Shuttle_FEs_SI(i,j) = Shuttle_FE;
        % Sat_fuels_SI(i,j) = Sat_FE;
        % leftoverShuttle_FMs_SI(i,j) = leftoverShuttle_FM;

        % Progress counter
        progress = int2str((i-1)*NE+j + 2*NA*NE) + "/" + int2str(3*NA*NE) + " analyses run"
    end
end



%% Plotting

% fakeNsats = 20;

% m_cutoff = 50; % kg
% Sat_masses(Sat_masses <= m_cutoff*fakeNsats) = NaN;

[X,Y] = meshgrid(eMOGs, aMOGs/AU);
figure(667)
surf(X,Y,Sat_masses)
% zscale log
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]");
zlabel("MOG Satellite Mass [kg]");
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", Hohmann Insertion")

figure(1667)
surf(X,Y,Sat_masses_1BI)
% zscale log
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]");
zlabel("MOG Satellite Mass [kg]");
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", 1-burn Insertion")

figure(2667)
surf(X,Y,Sat_masses_SI)
% zscale log
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]");
zlabel("MOG Satellite Mass [kg]");
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", Spread-out Insertion")

min_val = 1e2;
max_val = maxSSPayload2LEO/NSats;

figure(668)
% contourf(X, Y, Sat_masses/fakeNsats, 2000, LineColor='none')
x = [min(X) max(X)];
y = [min(Y) max(Y)];
h1 = imagesc(x,y,Sat_masses);
set(h1, 'AlphaData', 1-isnan(Sat_masses));
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", Hohmann Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("MOG eccentricity eMOG")
ylabel("aMOG [AU]")
ylim([0 6])

c = colorbar;
ylabel(c,'Dry mass per satellite [kg]', FontSize=12);
colormap turbo

figure(1668)
h2 = imagesc(x,y,Sat_masses_1BI);
set(h2, 'AlphaData', 1-isnan(Sat_masses_1BI));
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", 1-burn Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("MOG eccentricity eMOG")
ylabel("aMOG [AU]")
ylim([0 6])

c = colorbar;
ylabel(c,'Dry mass per satellite [kg]', FontSize=12);
colormap turbo

figure(1669)
h3 = imagesc(x,y,Sat_masses_SI);
set(h3, 'AlphaData', 1-isnan(Sat_masses_SI));
title("MOG Satellite Masses with NSats = " + num2str(NSats) + ", Spreadout Insertion")
set(gca,'ColorScale','log','YDir','normal', 'CLim', [min_val max_val]);
xlabel("MOG eccentricity eMOG")
ylabel("aMOG [AU]")
ylim([0 6])

c = colorbar;
ylabel(c,'Dry mass per satellite [kg]', FontSize=12);
colormap turbo

%% AAS/AIAA Figures

clf(669)
figure(669)
LEGEND = [];
hold on
n = 5;
for i = floor(linspace(NE/(n+1),NE*n/(n+1),n))
    plot(aMOGs/AU, Sat_masses(:,i))
    plot(aMOGs/AU, Sat_masses_1BI(:,i))
    LEGEND = [LEGEND, "eMOG = " + num2str(eMOGs(i)) + " (H)", "eMOG = " + num2str(eMOGs(i)) + " (1BI)"];
end
hold off
yscale log
legend(LEGEND)
xlabel("aMOG [AU]")
ylabel("MOG Satellite Mass [kg]")

figure(2832)
surf(X,Y,Sat_masses_1BI./Sat_masses)
title("Ratio of Mass per MOG Sat with NSats = " + num2str(NSats) + ", 1-burn Insertion/ Hohmann Insertion")
xlabel("eMOG")
ylabel("aMOG [AU]")
zlabel("Ratio")

figure(2833)
surf(X,Y,Sat_masses_SI./Sat_masses)
title("Ratio of Mass per MOG Sat with NSats = " + num2str(NSats) + ", Spreadout Insertion/ Hohmann Insertion")
xlabel("eMOG")
ylabel("aMOG [AU]")
zlabel("Ratio")

clf(28932)
figure(28932)
plot(min(Sat_masses_1BI./Sat_masses), 'red')
hold on
plot(max(Sat_masses_1BI./Sat_masses), 'red')
plot(min(Sat_masses_SI./Sat_masses), 'blue')
plot(max(Sat_masses_SI./Sat_masses), 'blue')
hold off
xlabel("eMOG")
ylabel("Ratio of Hohmann Transfer Satellite Mass to 1-Burn Insertion Satellite Mass")
legend("min 1BI/2BI", "max 1BI/2BI", "min SI/2BI", "max SI/2BI");





