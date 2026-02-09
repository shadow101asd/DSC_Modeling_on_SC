% PLotting an example FC for the Parametrization section

% Add paths

addpath("../../Run/") % Run folder
addpath("../") % Main Plotting scripts

% Constants

AU = 149600000; % km
muSu = 1.327124400419393e+11;

% Parameters

Nw = 2;
NM = 3;

aF = AU;
eF = 0.35;
w0F = pi/3;
f0F = pi/6;

Nsats = 17;
etR = 0;
res_outline = 500;

% Get FC

XSats = generate2DP_Flower(aF,eF,w0F,f0F,Nw,-NM,etR,muSu,Nsats);
XSats_outline = generate2DP_Flower(aF,eF,w0F,f0F,Nw,-NM,etR,muSu,res_outline);

% Plotting

t = 1;
lim = 1.15*aF*(1+eF)/AU;
n = 2359248;
f1 = figure(n);
clf(f1)
hold on
viscircles([0,0; 0,0], [(1-eF); (1+eF)]*aF/AU, LineStyle="-.", LineWidth=1, Color="blue")
viscircles([0,0], aF/AU, LineStyle="-", LineWidth=1, Color="blue")
h1 = scatter(0,0,SizeData=150, MarkerEdgeColor="black", MarkerFaceColor="yellow", Marker="pentagram", LineWidth=0.75);
h3 = plot(permute(XSats_outline(1,t,:)/AU, [3,1,2]),permute(XSats_outline(2,t,:)/AU, [3,1,2]), ...
    LineStyle="--", LineWidth=2, Color="black");
h2 = scatter(permute(XSats(1,t,:)/AU, [3,1,2]),permute(XSats(2,t,:)/AU, [3,1,2]), ...
    SizeData=120, MarkerEdgeColor="black", MarkerFaceColor=[0.7, 0.7, 0.7], Marker="square");
hold off
axis equal
xlim([-lim lim])
ylim([-lim lim])
xlabel("x [AU]")
ylabel("y [AU]")
legend([h1, h2, h3], "Sun", "FC Satellites", "FC Track")
set(gca, FontName="Times New Roman", FontSize=13)
set(f1, Position=[500 500 400 400])

%% Saving

filename = "./Figures/exampleFC_parametrization.eps";
exportgraphics(f1, filename, BackgroundColor='none')


