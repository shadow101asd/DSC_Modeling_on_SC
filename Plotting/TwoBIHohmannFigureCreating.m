
AU = 149600000;
muSu = 1.327124400419393e+11;
aMOG = AU;
eMOG = 0.5;
ata = aMOG*(1+eMOG/2);
atp = aMOG*(1-eMOG/2);
eta = (aMOG*eMOG)/(aMOG + aMOG*(1+eMOG));
etp = (aMOG*eMOG)/(aMOG + aMOG*(1-eMOG));

figure(124231)
clf(124231)
scatter(0, 0, "pentagram", "yellow", "filled")
hold on
plotEllipse([aMOG, 0, 0, 0, 0, 0], muSu, 1, AU, "black")
plotEllipse([aMOG, eMOG, 0, 0, 0, 0], muSu, 1, AU, "blue")
% Transfer Orbits
plotEllipse([ata, eta, 0, 0, 0, 0], muSu, 0.5, AU, "red")
plotEllipse([atp, etp, 0, 0, 0, 0], muSu, -0.5, AU, "green")
hold off
lim = 1.75;
xlim([-lim, lim])
ylim([-lim, lim])
xlabel("x [AU]")
ylabel("y [AU]")
legend("Sun", "Circular Orbit", "MOG Orbit", "Hohmann transfer to MOG Orbit Apoapsis", "Hohmann transfer to MOG Orbit Periapsis")