% MOG Velocities in Polar Non-inertial Frame Computation

etR = 1; % Doesn't matter here
AU = 1.496e8; % 1 AU in km
aMOG = 1.25*AU; % km
eMOG = 0.6;
muSu = 1.327124400419393e+11;
day = 3600*24; % a day in seconds

[theta, rho] = getPolarRepresentationOfMOG(aMOG, eMOG, muSu); % rho is in AU

rho = AU*rho; % Converting from AU to km

% get rid of discontinuity in theta
[~, idx] = max(theta);
theta(idx+1:end) = theta(idx+1:end) + 2*pi;

% Get array of corresponding time vectors
TMOG = 2*pi*sqrt(aMOG^3/muSu);
ts = linspace(0, TMOG, length(rho)+1)';

% Pad end of vectors for velocity computations (periodic behavior)
theta(end+1) = theta(1);
rho(end+1) = rho(1);

% Compute velocities along radial and tangential directions

angularVs = mod((theta(2:end) - theta(1:end-1)), 2*pi)./(ts(2:end) - ts(1:end-1)); % radians/second
radialVs = (rho(2:end) - rho(1:end-1))./(ts(2:end) - ts(1:end-1)); % km/s

%% Fitting

% % Rho
% [fit_rho, gof_rho] = fit(ts(1:end-1)/day, rho(1:end-1)/AU, "fourier8");
% 
% % Theta
% theta_lin = linspace(0,2*pi,length(ts(1:end-1)))';
% [fit_theta, gof_theta] = fit(ts(1:end-1)/day, theta(1:end-1)-theta_lin, "fourier8");

fourierDegree = 8;
resolution = 1e4;
[fit_rho, fit_theta, gof_rho, gof_theta, fake_t] = numericMOGShapeFitting(eMOG, fourierDegree, resolution);

% Analysis

[xreal, yreal] = pol2cart(theta(1:end-1), rho(1:end-1)); % km and rad
[xfit, yfit] = pol2cart(fit_theta(fake_t), fit_rho(fake_t) * eMOG*aMOG); % km and rad

RMSE = sqrt(mean((xreal - xfit).^2) + mean((yreal - yfit).^2));

% Calculate Error
E = sqrt((xreal-xfit).^2 + (yreal-yfit).^2);

MAE = mean(E);

MAE_percentaMOG = MAE/aMOG *100 % [%]

%% Plotting

figure(998)
plot(ts(1:end-1)/day, rho(1:end-1)/AU)
hold on
plot(ts(1:end-1)/day, fit_rho(fake_t) * eMOG*aMOG)
hold off
xlabel("time [days]");
ylabel("distance from center of MOG [AU]")
legend("Data", "Fourier fit");

figure(999)
plot(ts(1:end-1)/day, theta(1:end-1))
hold on
plot(ts(1:end-1)/day, fit_theta(fake_t))
hold off
xlabel("time [days]");
ylabel("angular position over MOG [rad]")
legend("Data", "Fourier fit");

figure(1000)
plot(ts(1:end-1)/day, radialVs)
xlabel("time [days]");
ylabel("radial velocity [km/s]")

figure(1001)
plot(ts(1:end-1)/day, angularVs*day*1000)
xlabel("time [days]");
ylabel("angular velocity [mrad/day]")

figure(1002)
polarplot(theta, rho/AU)
hold on
polarplot(fit_theta(fake_t), fit_rho(fake_t) * eMOG*aMOG/AU)
hold off
legend("Data", "Composite Fourier fit");

figure(93849357)
plot(ts(1:end-1)/day, xreal/AU);
hold on
plot(ts(1:end-1)/day, yreal/AU)
plot(ts(1:end-1)/day, xfit/AU)
plot(ts(1:end-1)/day, yfit/AU)
hold off
xlabel("time [days]");
ylabel("[AU]")
legend("x data", "y data", "x estimate", "y estimate");

figure(100)
plot(ts(1:end-1)/day, E);
xlabel("time [days]");
ylabel("error [km]");
title("Error in Cartesian Coordinates");

