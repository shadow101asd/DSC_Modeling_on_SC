% PLotting 2D Flower Constellation Phasing DV

AU = 149600000;
muSu = 1.327124400419393e+11;

W = 3;
M = -2;
eF = 0.4;

resolution = 100;
dPs = linspace(0,2*pi,resolution);

Nm = 1;
[DVs, DVs_seq] = computePhasingDVs_2DF_sequentialM(W, M, eF, resolution, Nm);

%% Dedicated Deployment Sim

aF = AU;
L = 5;
Npl = 5;
frac = 1;
[Sat_DVs, DPs_opt, pattern_opt, maxDV] = compute_SIdeployment_SatDVs_2DF(W,M,aF,eF,L,Npl,frac,resolution);

patterns = slots_1D_tile_nonperiodic_recursive(L,Npl);
[DVs_of_pattern_curve, DPs_for_pattern_curve] = get_DVs_for_pattern(DVs, patterns, frac);


%% Plotting

plot_resolution = 500;
XSats = generate2DP_Flower(AU,eF,pi/2,0,W,M,1,muSu,plot_resolution);
limval = ceil(2*max(max(max(abs(XSats))))/AU + 0.01)/2; % AU
figure(123)
scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
hold on
scatter(reshape(XSats(1,1,:), [size(XSats, 3),1])/AU,reshape(XSats(2,1,:), [size(XSats, 3),1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
xlabel("x [AU]");
ylabel("y [AU]");
hold off
title("W = " + W + " | M = " + M + " | N_{sats} = " + plot_resolution)
set(gca, FontSize=14)

figure(1213)
plot(dPs,DVs, LineWidth=2)
hold on
plot(dPs,DVs_seq, LineWidth=2)
plot(DPs_for_pattern_curve,DVs_of_pattern_curve, LineWidth=2)r
scatter(DPs_opt, Sat_DVs, Marker="+", SizeData=100, MarkerEdgeColor="black", LineWidth=2)
hold off
xlabel("Flower Constellation Phase Change (rad)")
ylabel("Phasing DV (km/s)")
set(gca, FontSize=14)