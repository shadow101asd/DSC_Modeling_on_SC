% Plotting the example of multiple eccentricities figure for the MOG
% presentation


AU = 149600000; % 1 AU in km
outline_res = 500;
load("../Inputs/run010.mat", "muSu");


% Define MOG parameters to be represented

aMOG = 1.0*AU;
eMOGs = [0.1, 0.4, 0.7];
w0s = [pi/8, pi/2, -4.2*pi/6];

% Generate Ephemerides

K1 = [aMOG 0 0 0 0 w0s(1)];
K2 = [aMOG 0 0 0 0 w0s(2)];
K3 = [aMOG 0 0 0 0 w0s(3)];

X1 = generateMOGNearPlanet(Keplerian2Cartesian(K1, muSu), muSu, 1, eMOGs(1), 1, outline_res-1, false);
X2 = generateMOGNearPlanet(Keplerian2Cartesian(K2, muSu), muSu, 1, eMOGs(2), 1, outline_res-1, false);
X3 = generateMOGNearPlanet(Keplerian2Cartesian(K3, muSu), muSu, 1, eMOGs(3), 1, outline_res-1, false);

% Pad to close the outline
X1(:,1,outline_res) = X1(:,1,1);
X2(:,1,outline_res) = X2(:,1,1);
X3(:,1,outline_res) = X3(:,1,1);

%% Plotting

figure(133234234)

scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);

limval = 1.8;
set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);

hold on
viscircles([0,0], 1, Color="blue", LineWidth=1.5);
plot(12, 12, Color="blue", LineWidth=1.5);
plot(reshape(X1(1,1,:), [outline_res,1])/AU,reshape(X1(2,1,:), [outline_res,1])/AU, LineWidth=3, Color="red");
plot(reshape(X2(1,1,:), [outline_res,1])/AU,reshape(X2(2,1,:), [outline_res,1])/AU, LineWidth=3, Color="green");
plot(reshape(X3(1,1,:), [outline_res,1])/AU,reshape(X3(2,1,:), [outline_res,1])/AU, LineWidth=3, Color="black");
hold off

xlabel('x [AU]');
ylabel('y [AU]');
legend("Sun", "Circle of radius 1 AU", "e_{MOG} = " + eMOGs(1), "e_{MOG} = " + eMOGs(2), "e_{MOG} = " + eMOGs(3))