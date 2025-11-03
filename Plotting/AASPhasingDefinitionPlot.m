

AU = 149600000;
muSu = 1.327124400419393e+11;
eMOG = 0.6;
aMOG = AU;
resolution = 5000;
cross_resolution = 8;

X1s = generate2DP_Flower(AU, eMOG, pi/2, 0, 1, -1, 1, muSu, resolution);
X1s = reshape(X1s, [6,resolution]);

X2s = generate2DP_Flower(AU, eMOG, pi/2, 0, 1, -1, 1, muSu, cross_resolution);
X2s = reshape(X2s, [6,cross_resolution]);

% Plotting

figure(23231)
plot(X1s(1,:)/AU, X1s(2,:)/AU, LineWidth=2, Color="black")
hold on
scatter(0,1, MarkerEdgeColor="red", Marker="+", LineWidth=2, SizeData=100)
scatter(X2s(1,:)/AU, X2s(2,:)/AU, MarkerEdgeColor="black", Marker="+", LineWidth=2.5, SizeData=100)
hold off
xlim([-1.25 1.25])
ylim([0 1.8])
pbaspect([2.5 1.8 1])
axis off


%% Polar plot

[theta1s, rho1s] = cart2pol(X1s(1,:), X1s(2,:)-aMOG);
[theta2s, rho2s] = cart2pol(X2s(1,:), X2s(2,:)-aMOG);

% normalize rhos
rho1s = rho1s / (aMOG*eMOG);
rho2s = rho2s / (aMOG*eMOG);

figure(23232)
polarplot(theta1s,rho1s,LineWidth=2, Color="black");
hold on
polarscatter(0,0,MarkerEdgeColor="red", Marker="+", LineWidth=2, SizeData=100)
polarscatter(theta2s,rho2s,MarkerEdgeColor="black", Marker="+", LineWidth=2.5, SizeData=100)
hold off