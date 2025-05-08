%% Numerical Determination of MOG Relative Shape

etR = 1;
AU = 1.496e8; % 1 AU in km
aMOG = 1.25*AU;
mu = 1; % Value doesn't matter

% eMOGs = linspace(1e-3, 1-1e-3, 5);
eMOG = 0.8;

% clf(50)
% figure(50)
% polarplot(0, 0);
% hold on
% for i = 1:length(eMOGs)
%     eMOG = eMOGs(i);
%     [theta, rho] = getPolarRepresentationOfMOG(aMOG/(eMOG), eMOG, mu);
%     polarplot(theta, rho);
% end
% hold off

% % Custom fit type
% ft = fittype('sqrt(am/(1-bm*cos(x)^2))*(a1*cos(x*p)+b1*sin(x*p)+a2*cos(2*x*p)+b2*sin(2*x*p))');
% coeffnames(ft)
% options = fitoptions(ft);
% options.StartPoint = [1, 1, 1, 1, 1,sqrt(3)/2, 1];
% options.Lower = [-1, -1, 1, -1, -1, sqrt(3)/2, -1];
% options.Upper = [ 1,  1, 1,  1,  1, sqrt(3)/2,  1];

% [theta, rho] = getPolarRepresentationOfMOG(aMOG/(eMOG), eMOG, mu);
% f1 = fit(theta,rho,"fourier8");
% N_samples = length(theta);
% % theta_even = linspace(0, 2*pi, N_samples);
% theta_even = theta;

[theta, rho] = getPolarRepresentationOfMOG(aMOG, eMOG, mu);

figure(51)
polarplot(theta, rho/min(rho), LineWidth=3);
% hold on
% polarplot(theta_even, f1(theta_even))
% hold off
% legend("Data", "Fourier8 Fit");