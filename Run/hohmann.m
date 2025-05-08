% Hohmann plot

AU = 1.496e8; % 1 AU in km
R1 = AU;
muSu = 1.327124400419393e+11;

R2s = AU*linspace(0.5, 1.6, 200);
SS_DVs = NaN(size(R2s));

for i = 1:length(R2s)
    R2 = R2s(i);
    SS_DVs(i) = computeHohmannTransferCirc2Circ(R1, R2, muSu);
end


%% Plotting

figure(100)
plot(R2s/AU, abs(SS_DVs))
xlabel("R2 [AU]");
ylabel("Magnitude of total DV for Hohmann Transfer [km/s]")