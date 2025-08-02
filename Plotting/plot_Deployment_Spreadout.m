%% Spread-out Insertion Plotting for AAS/AIAA MOG Presentation

% Parameters

AU = 149600000; % 1 AU in km

outline_res = 500;

interval = 3600*24*1; % 1 day
sim_end = 3600*24*365.25*4.5; % 4.5 years
etR = 0:interval:sim_end; % Row vector of times between start and end date

load("../Inputs/run010.mat", "muSu")

KEa = [AU; 0; 0; 0; 0; 0];
XEa = propagateFromKeplerians(KEa, muSu, etR);

% Define MOG parameters to be represented

aMOG = 1.4*AU;
eMOG = 0.25;
Nsats = 7;

% Representation indices

shuttle_launch_idx = floor(length(etR)/16);

%% Generate shuttle trajectory

XShuttle_E = XEa(:,1:shuttle_launch_idx);

% Shuttle transfer orbit (circ 2 MOG ell here)

at = 0.5*(AU+aMOG*(1+eMOG));
f0 = 0;
et = (aMOG*(1+eMOG)-AU)/(aMOG*(1+eMOG)+AU);
[~,~,~,~,~,wt] = Cartesian2Keplerian(XEa(:,shuttle_launch_idx), muSu);

K_st_i = [at, et, 0, 0, wt, f0];
t_transfer = pi*sqrt(at^3 / muSu);
shuttle_arrival_idx = shuttle_launch_idx+1+ floor(t_transfer/interval);
etR_transfer = etR(shuttle_launch_idx+1:shuttle_arrival_idx);

XShuttle_transfer = propagateFromKeplerians(K_st_i, muSu, etR_transfer);

% Shuttle deployment orbit

etR_deploy = etR(shuttle_arrival_idx+1:end);
K_sd_i = [aMOG; eMOG; 0; 0; wt; pi];
XShuttle_deploy = propagateFromKeplerians(K_sd_i, muSu, etR_deploy);

XShu = cat(2, XShuttle_E, XShuttle_transfer, XShuttle_deploy);

%% Generate MOG Satellite Trajectories

XSats = NaN([6, length(etR), Nsats]); % Sats default to not appearing

for n = 1:Nsats
    % Transfer
    
end

end_of_deployment_idx = transfer_end_idx;

%% Generate MOG Outline for Plotting
% 
% dtheta = pi*((1 + eMOG/2)^(3/2) - 1); % see paper
% KO = [aMOG, 0, 0, 0, 0, (K_sd_i(6)-dtheta)];
% XO = propagateFromKeplerians(KO, muSu, [0, -etR(shuttle_arrival_idx)]);
% XO = XO(:,2);
% XOutline = generateMOGNearPlanet(XO, muSu, etR, eMOG, 1, outline_res, false);

%% Setup video export

% Create a high-quality MP4 video from animated plot

% v = VideoWriter('PlottingData/AAS Paper Figs and Data/MOG_Deployment_SI_pres.mp4', 'MPEG-4');
% v.Quality = 100;        % Max quality (0–100)
% v.FrameRate = 90;       % Adjust to match your desired playback speed
% open(v);

%% Plotting

fig_num = 3242344;
fig = figure(fig_num);

limval = 2; % AU
for t = 1:length(etR)
    clf(fig_num)

    scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black');
    set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
    
    hold on
    % viscircles([0,0], aMOG/AU, Color="black", LineWidth=1); % TODO:
    % replace with Shuttle deployment orbit?
    viscircles([0,0], 1, Color="blue", LineWidth=1);

    scatter(XEa(1,t)/AU,XEa(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="blue", Marker="o", SizeData=350);
    scatter(XShu(1,t)/AU,XShu(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="black", Marker="^", SizeData=350);
    scatter(reshape(XSats(1,t,:), [Nsats,1])/AU,reshape(XSats(2,t,:), [Nsats,1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
    % plot(reshape(XOutline(1,t,:)/AU, [outline_res, 1]), reshape(XOutline(2,t,:)/AU, [outline_res, 1]), LineStyle="--", LineWidth=2, Color='black')

    % Add axis labels and title
    xlabel('x [AU]');
    ylabel('y [AU]');
    
    if t < end_of_deployment_idx
        title("Time since start of deployment procedure: " + num2str((etR(t)-etR_transfer(1))/(3600*24*365.25), '%.2f') + " Earth years")
    else

        title("Total deployment time: " + num2str((etR(end_of_deployment_idx)-etR_transfer(1))/(3600*24*365.25), '%.2f') + " Earth years")
    end
    
    hold off
    pause(1e-2);
    % exportgraphics(gcf,'PlottingData/AAS Paper Figs and Data/MOG_Deployment_Hohmann_pres.gif','Append',true);
    
    % Export as image (rasterized, clean rendering)
    % exportgraphics(fig, 'temp_frame.png', 'Resolution', 300);
    % frame = imread('temp_frame.png');
    % writeVideo(v, frame);
end
% 
% close(v)
% delete('temp_frame.png');