%% Running
tic

dPs = linspace(0.0, 2*pi, 1000);

DVs = Inf(size(dPs));
ToPs = DVs;
muSu = 1.327124400419393e+11;
AU = 149600000; % km

aMOG = 1.0*AU;
eMOG = 0.1;

for i = 1:length(dPs)
    dP = dPs(i);
    [DV, ToP] = computeOptimalPhasingDVMOG(aMOG, eMOG, dP, muSu);
    DVs(i) = DV;
    ToPs(i) = ToP;
    disp("Progress: " + num2str(i) + "/" + num2str(length(dPs)) + " analyses")
end

toc

%% Analysis
% options = fitoptions( 'Method', 'NonlinearLeastSquares',...
%     'Lower',[0], 'Upper',[100], 'StartPoint', 5);
% fit1 = fit(dPs(2:end)', DVs(2:end)', 'a*x')

fit1 = fit(dPs(1:end)', DVs(1:end)', 'poly2')

%% Plotting

figure(6666)
plot(dPs, DVs)
hold on
% plot(dPs, fit1(dPs))
hold off
xlabel("MOG Phase Shift [rad]");
ylabel("Optimal DV Found [km/s]");
% legend("Data", "Fit")
title("")

day = 24*3600;
figure(6667)
plot(dPs, ToPs/day)
xlabel("MOG Phase Shift [rad]");
ylabel("Time of Flight for Phase Shift (Earth days)");

figure(6668)
plot(dPs, DVs, LineWidth=2)
hold on
scatter(dPs(1,[1,333,667]), DVs(1,[1,333,667]), 'filled', SizeData=200, Marker='diamond')
xline(dPs(1,[1,333,667]), LineStyle="--", Color="black", LineWidth=1.5);
hold off
xlim([0, 2*pi])
xlabel("MOG Phase Shift Δϕmog [rad]");
ylabel("Optimal ΔV Found [km/s]");
legend("ΔV", 'Selected ΔVi', "Selected Δϕmog")
% legend("Data", "Fit")

figure(6669)
plot(dPs, DVs, LineWidth=2)
hold on
scatter(dPs(1,[125,375,625,875]), DVs(1,[125,375,625,875]), 'filled', SizeData=200, Marker='diamond')
xline(dPs(1,[125,375,625,875]), LineStyle="--", Color="black", LineWidth=1.5);
hold off
xlim([0, 2*pi])
xlabel("MOG Phase Shift Δϕmog [rad]");
ylabel("Optimal ΔV Found [km/s]");
legend("ΔV", 'Selected ΔVi', "Selected Δϕmog")
% legend("Data", "Fit")
