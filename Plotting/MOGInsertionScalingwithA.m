% Evaluating the Scaling w.r.t. aMOG of MOG insertion maneuvers

AU = 149600000;
muSu = 1.327124400419393e+11;
day = 3600*24;

%% Run

eMOG = 0.3;
Na = 200;
aMOGs = linspace(0.2,1.8,Na)'*AU;

DVs = zeros(size(aMOGs));
TOFs = zeros(size(aMOGs));

for i = 1:Na
    [MOGSat_DV, LToF] = computeHohmannTransferCirc2Ell(aMOGs(i), aMOGs(i), eMOG, muSu);
    DVs(i) = MOGSat_DV;
    TOFs(i) = LToF;
end

%% Fitting

[fitDV, gofDV] = fit(aMOGs/AU, DVs, "a * x^(-1/2)") % This follows from the aMOG scaling as expected
[fitTOF, gofTOF] = fit(aMOGs/AU, TOFs/day, "a * x^(3/2)") % This one is just half of the period of the Hohmann transfer orbit


%% Plotting

figure(1)
plot(aMOGs/AU, DVs)
xlabel("aMOG [AU]")
ylabel("MOG Satellite Insertion DV [km/s]")

figure(20)
plot(aMOGs/AU, TOFs/day)
xlabel("aMOG [AU]")
ylabel("MOG Satellite Insertion Time of Flight [days]")
