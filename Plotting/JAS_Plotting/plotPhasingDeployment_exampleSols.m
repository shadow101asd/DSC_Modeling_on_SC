% Plotting the FC Phase Change Example Results

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

aF = 1*AU;
eF = 0.5;

Nw = 3;
NM = -2;

NL = 3;
Npl = 4;
frac = 0.5;

resolution = 1000;

%% Run

patterns = slots_1D_tile_nonperiodic_recursive(NL,Npl); % Possible DP patterns from combinatorics problem

DVs = computePhasingDVs_2DF_sequentialM(Nw, -NM, eF, resolution, 1)*sqrt(AU/aF); % scale correctly according to aF
[DVs_pats, DPs_pats] = get_DVs_for_pattern(DVs, patterns, frac);

%% Find best sol (manual stuff for sure here)

[minDV, xopt] = min(DVs_pats(2,:));

offsets = linspace(0,2*pi,length(DVs_pats));
dphi_opts = mod(2*pi*frac*[0, 1, 6, 7]/12 + offsets(xopt), 2*pi);
DVs_opts = zeros(1, Npl);
for i = 1:Npl
    DVs_opts(i) = computeOptimalPhasingDV_2DF_extended(Nw, NM, eF, dphi_opts(i));
end

%% Plotting

n = 987656343;
w = 600;
h = 250;

f1 = figure(n);
plot(linspace(0,2*pi,length(DVs)), DVs', LineWidth=1)
hold on
scatter(dphi_opts, DVs_opts, "+", LineWidth=2, SizeData=100, MarkerEdgeColor="black")
hold off
xlim([0, 2*pi])
ylim([0 25])
xlabel("\Delta\phi_F [radians]")
ylabel("\DeltaV [km/s]")
set(gca, FontName="Times New Roman", FontSize=13)
set(f1, Position=[300 500 w h])

f2 = figure(n+1);
plot(offsets, DVs_pats', LineWidth=1)
hold on
scatter(offsets(xopt), minDV, "+", LineWidth=2, SizeData=100, MarkerEdgeColor="black")
hold off
xlim([0, 2*pi])
ylim([0 25])
set(f2, Position=[300+w 500 w h])
legend("P = [1 0 0 1 0 0 1 0 0 1 0 0]", ...
       "P = [1 1 0 0 0 0 1 1 0 0 0 0]", ...
       "P = [1 1 1 1 0 0 0 0 0 0 0 0]", ...
       Location="best")
xlabel("\Delta\phi_F^1 [radians]")
ylabel("\DeltaV [km/s]")
set(gca, FontName="Times New Roman", FontSize=13)

%% Saving

filename1 = "./Figures/PhasingDeployment_Methodology_Plot1.eps";
exportgraphics(f1, filename1, BackgroundColor="white")

filename2 = "./Figures/PhasingDeployment_Methodology_Plot2.eps";
exportgraphics(f2, filename2, BackgroundColor="white")


