% Plotting the DD SatDV example curves

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

aF = 1*AU;

E = 100;
eFs = linspace(1e-3, 0.8, E);

Nws = [1, 3, 1, 21];
NMs = [1, 2, 2, 6];
A = length(NMs);

%% Run

DVs = NaN([A,E]);
resolution = 10000;

for i = 1:E
    for a = 1:A
        DVs(a,i) = compute_DIdeployment_SatDV_2DF(Nws(a),NMs(a),aF,eFs(i),resolution);
    end
end


%% Plotting

n = 33876;

labels = [];
for a = 1:A
    labels = [labels, "(N_\omega, N_M) = (" + int2str(Nws(a)) + "," + int2str(NMs(a)) + ")"];
end

colors = [
    0.00 0.45 0.70
    0.90 0.60 0.00
    0.80 0.60 0.70
    0.00 0.60 0.50
    0.80 0.40 0.00
    0.00 0.00 0.00];
patterns = ["--", "-.", "-", ":"];
markers = ["o", "s", "^", "d"];

w = 500;
h = 300;

f1 = figure(n);
clf(f1)
hold on
for a = 1:A
    plot(eFs, DVs(a,:), Color=colors(a,:), LineStyle="-", Marker=markers(a), ...
        MarkerIndices=(1:floor(E/6):E), MarkerFaceColor='w', LineWidth=1)
end
hold off
legend(labels, Location="best")
xlabel("FC eccentricity e_F")
ylabel("\DeltaV [km/s]")
set(gca, FontName="Times New Roman", FontSize=13)
set(f1, Position=[200 200 w h])


%% Saving

filename = "./Figures/DirectDeployment_exampleSatDVs.eps";
exportgraphics(f1, filename, BackgroundColor="white")