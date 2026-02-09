% Plotting the FC Phase Change Example Results

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

eFs = [0.1, 0.3, 0.5];
E = length(eFs);

Nws = [1 3 2  1];
NMs = [1 2 5 -1];
assert(length(NMs) == length(Nws))
F = length(Nws);

resolution = 500;
eps = 1e-2;
dPs = linspace(eps, 2*pi-eps, resolution);

% Initialize fields
DVs = zeros(E,resolution,F);

for f = 1:F
    Nw_f = Nws(f);
    NM_f = NMs(f);
    for e = 1:E
        eF_e = eFs(e);
        parfor i = 1:resolution
            DVs(e,i,f) = computeOptimalPhasingDV_2DF_extended(Nw_f, -NM_f, eF_e, dPs(i));
        end
    end
end


%% Plotting

n = 7752448;
h = 300; % pixels
w = 500; % pixels
screen_h = 700;
screen_w = 1500; % pixels

colors = [
    0.00 0.45 0.70
    0.90 0.60 0.00
    0.80 0.60 0.70
    0.00 0.60 0.50
    0.80 0.40 0.00
    0.00 0.00 0.00];
patterns = ["--", "-.", "-", ":"];
markers = ["o", "s", "^", "d"];

fs = [];
for f = 1:F
    fs(f) = figure(n+f);
    clf(fs(f))
    hold on
    for e = 1:E
        plot(dPs, DVs(e,:,f), LineWidth=1, Color=colors(e,:), LineStyle="-", Marker=markers(e), MarkerIndices=(16:50:resolution), MarkerFaceColor='w');
    end
    hold off
    xlabel("\Delta\phi_F [radians]")
    ylabel("\DeltaV [km/s]")
    legend(arrayfun(@(x) sprintf('e_F = %.3g', x), eFs, 'UniformOutput', false), Location="best");
    xlim([0 2*pi])              % adjust as needed
    xticks(0:pi/2:2*pi)
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
    set(gca, FontName="Times New Roman", FontSize=13)
    set(fs(f), Position=[(mod(w*(f-1),screen_w)) (screen_h-h*floor(w*(f-1)/screen_w)) w h])
end

%% Saving

for i = 1:F
    filename = "./Figures/Phasing_results_" + int2str(i) + ".eps";
    exportgraphics(figure(fs(i)), filename, BackgroundColor="white")
end