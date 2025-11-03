%% Flower Constellation Phasing Maneuver Analysis

AU = 149600000;
muSu = 1.327124400419393e+11;

W = 3;
M = -2;
Npl = 5;
frac = 1;

NeF = 25;
eFs = linspace(1e-2, 0.5, NeF);

NL = 25;
Ls = 1:NL;

worstCaseSatDVs = NaN([NL,NeF]);
tic
parfor L = Ls
    for j = 1:NeF
        eF = eFs(j);
        [~,~,~,maxDV] = compute_SIdeployment_SatDVs_2DF(W,M,aF,eF,L,Npl,frac,resolution);
        worstCaseSatDVs(L,j) = maxDV;
    end
end
toc

%% Plotting

[X,Y] = meshgrid(eFs, Ls);
x = [min(X) max(X)];
y = [min(Y) max(Y)];

figure(121243)
h1 = imagesc(x,y,worstCaseSatDVs);
title("Worst-case Satellite DV of F2D Spread-out Insertion Maneuver, W = " + W + ", M = " + M + ", Npl = " + Npl)
set(gca,'YDir','normal');
xlabel("Flower Constellation eccentricity eF")
ylabel("Number of launches to deploy the Flower Constellation")

c = colorbar;
ylabel(c,'Worst-case Satellite DV [km/s]', FontSize=12);
colormap parula