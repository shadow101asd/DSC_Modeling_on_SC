% Hohmann DV Comparison Plotting for AAS/AIAA MOG Paper


AU = 149600000;
aMOG = AU;
muSu = 1.327124400419393e+11;
epsilon = 1e-4;
eMOGs = linspace(epsilon, 1-epsilon, 1/epsilon-1);

DV_as = NaN(size(eMOGs));
DV_ps = NaN(size(eMOGs));

%% Run

for i = 1:length(eMOGs)
    [DV, ~, ~, ~, ~, worseDV] = computeHohmannTransferCirc2Ell(aMOG,aMOG,eMOGs(i),muSu);
    DV_as(i) = DV;
    DV_ps(i) = worseDV;
end

%% Plot

figure(2048320)
plot(eMOGs, DV_as)
hold on
plot(eMOGs, DV_ps)
plot(eMOGs, DV_ps-DV_as)
hold off
xlabel("MOG eccentricity, eMOG")
ylabel("DV, km/s")
legend("DVa", "DVp", "DVp - DVa")