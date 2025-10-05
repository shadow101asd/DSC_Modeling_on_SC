% Plot Feature Examples for IEEE Paper

% Parameters
AU = 149600000; % 1 AU in km
NSats = 50;
etR = 1;
lim = 1.75;

% Load XEa and XMa
run_idx = '017';
input_file = "../Inputs/run017.mat";
load(input_file, "muSu", "XEa", "XMa")
t = 1;
XEa = XEa(:,t);
XMa = XMa(:,t);

aMOG = 1.25*AU;
eMOG = 0.35;
w0 = -0.7;

% MOG
Ki = [aMOG; 0.25; 0.0; 0.0; w0; 0.0]; % Second element (eccentricity) isn't actually used here
XMOG = NSATSpropagateFromKepleriansSHELLS(Ki,muSu,etR,NSats,eMOG,1);

fignum1 = 3252452;
figure(fignum1)
plot_DSC_Modular_Structure(XMOG, run_idx, fignum1, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')

% MOG @ Planet
eMOG2 = 0.24;
XMOGatP = generateMOGNearPlanet(XMa, muSu, etR, eMOG2, 1, NSats, true);

fignum2 = fignum1+1;
figure(fignum2)
plot_DSC_Modular_Structure(XMOGatP, run_idx, fignum2, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')


% MOG including Planet
XMOGincP = generateMOGIncludingPlanet(XMa, muSu, etR, NSats);
fignum3 = fignum2+1;
figure(fignum3)
plot_DSC_Modular_Structure(XMOGincP, run_idx, fignum3, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')


% Circular Ring
aC = 1.2*AU;
f0C = 2.3;
fracC = 0.7;
XCirc = NSATSpropagateFromKepleriansCIRC(aC,f0C,muSu,etR,NSats,fracC);

fignum4 = fignum3+1;
figure(fignum4)
plot_DSC_Modular_Structure(XCirc, run_idx, fignum4, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')

% Elliptical Ring
aE = 1.0*AU;
eE = 0.6;
wE = 1.3;
f0E = 0;
fracE = 0.7;

Xcenter = Keplerian2Cartesian(aE,eE,0,0,wE,f0E,muSu);
XEll = generateXsAlongOrbit(Xcenter, muSu, etR, NSats, fracE);

fignum5 = fignum4+1;
figure(fignum5)
plot_DSC_Modular_Structure(XEll, run_idx, fignum5, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')


% Ring Including Earth

fracRIncP = 0.7;
XRIncP = generateXsAlongOrbit(XEa, muSu, etR, NSats, fracRIncP);

fignum6 = fignum5+1;
figure(fignum6)
plot_DSC_Modular_Structure(XRIncP, run_idx, fignum6, 'Mars');
xlim([-lim lim])
ylim([-lim lim])
set(gca,'XColor', 'none','YColor','none')

