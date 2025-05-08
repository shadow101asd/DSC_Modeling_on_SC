% MOG 1-Burn Insertion Intersection Plotting

AU = 149600000; % km
muSu = 1.327124400419393e+11; 
aMOG = 1.0*AU;
eMOG = 0.5;

figure(1297392)
clf(1297392)
scatter(0, 0, "pentagram", "yellow", "filled")
hold on
plotEllipse([aMOG, 0, 0, 0, 0, 0], muSu, 1, AU, "black"); % Circular Orbit
plotEllipse([aMOG, eMOG, 0, 0, -pi, 0], muSu, 1, AU, "blue"); % MOG Eliptical Orbit
% Intersections
theta = atan(sqrt(1-eMOG^2)/eMOG);
scatter([aMOG*cos(theta)/AU, aMOG*cos(-theta)/AU], [aMOG*sin(theta)/AU, aMOG*sin(-theta)/AU])
hold off
lim = 1.5;
xlim([-lim lim])
ylim([-lim lim])
xlabel("x [AU]")
ylabel("y [AU]")
legend("Sun", "Circular Orbit", "MOG Elliptical Orbit", "Intersections");
pbaspect([1 1 1])