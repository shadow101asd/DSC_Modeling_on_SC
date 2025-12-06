%% 1-Burn Insertion Plotting for AAS/AIAA MOG Presentation

% Parameters

AU = 149600000; % 1 AU in km

outline_res = 500;

interval = 3600*24*0.5; % 1 day
sim_end = 3600*24*365.25*5.0; % 3.2 years
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

% Shuttle transfer orbit (circ 2 circ here)

at = 0.5*(AU+aMOG);
f0 = 0;
et = (aMOG-AU)/(aMOG+AU);
[~,~,~,~,~,wt] = Cartesian2Keplerian(XEa(:,shuttle_launch_idx), muSu);

K_st_i = [at, et, 0, 0, wt, f0];
t_transfer = pi*sqrt(at^3 / muSu);
shuttle_arrival_idx = shuttle_launch_idx+1+ floor(t_transfer/interval);
etR_transfer = etR(shuttle_launch_idx+1:shuttle_arrival_idx);

XShuttle_transfer = propagateFromKeplerians(K_st_i, muSu, etR_transfer);

% Shuttle deployment orbit

etR_deploy = etR(shuttle_arrival_idx+1:end);
K_sd_i = [aMOG; 0; 0; 0; 0; wt+pi];
XShuttle_deploy = propagateFromKeplerians(K_sd_i, muSu, etR_deploy);

XShu = cat(2, XShuttle_E, XShuttle_transfer, XShuttle_deploy);

%% Generate MOG Satellite Trajectories

% XSats_N = NaN([6, shuttle_arrival_idx, Nsats]); % Sats don't appear during these phases
XSats_N = repmat(XShu(:,1:shuttle_arrival_idx), [1,1,Nsats]); % Sats default to being (visually) in the shuttle


% MOG
XS = [XShuttle_deploy(1:2,1); 1e-2; XShuttle_deploy(4:5,1); 1e-2];
XSats_MOG = generateMOGNearPlanet(XS, muSu, etR_deploy, eMOG, 1, Nsats, true);

% Hide Satellites until they're deployed

Ds = distanceBetweenXs(XShuttle_deploy, XSats_MOG);

appear_idxs = ones([Nsats, 1]);
for i = 1:Nsats
    appear_idxs(i) = find(islocalmin(Ds(:,:,i)), 1, "first");
    % XSats_MOG(:,1:appear_idxs(i), i) = NaN;
    XSats_MOG(:,1:appear_idxs(i), i) = XShu(:,(shuttle_arrival_idx+1):shuttle_arrival_idx+appear_idxs(i));
end

end_of_deployment_idx = max(appear_idxs) + shuttle_arrival_idx;

% Recombine
XSats = cat(2, XSats_N, XSats_MOG);

%% Generate MOG Outline for Plotting

XO = propagateFromKeplerians(K_sd_i, muSu, [0, -etR(shuttle_arrival_idx)]);
XO = XO(:,2);
XOutline = generateMOGNearPlanet(XO, muSu, etR, eMOG, 1, outline_res, true);

%% Plotting

fig_num = 3242343;
fig = figure(fig_num);
set(fig, 'Visible', 'off', 'Units', 'pixels', 'Position', [200 200 700 600], 'Resize', 'off');
% avoid GUI rendering issues

% Status indicators
stage(length(etR)) = "MOG Delivered";
stage(1:shuttle_launch_idx) = "Awaiting Launch";
stage((shuttle_launch_idx+1):shuttle_arrival_idx) = "Shuttle Transfering";
stage((shuttle_arrival_idx+1):end_of_deployment_idx) = "Satellites Deploying";
stage(end_of_deployment_idx:end) = "MOG Deployed";

background(length(etR)) = "MOG Delivered";
background(1:shuttle_launch_idx) = "blue";
background((shuttle_launch_idx+1):shuttle_arrival_idx) = "black";
background((shuttle_arrival_idx+1):end_of_deployment_idx) = "#AAAAAA";
background(end_of_deployment_idx:end) = "#339933";

limval = 2; % AU
for t = 1:length(etR)
    clf(fig_num)

    scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
    
    hold on
    viscircles([0,0], aMOG/AU, Color="black", LineWidth=1);
    viscircles([0,0], 1, Color="blue", LineWidth=1);

    scatter(XEa(1,t)/AU,XEa(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="blue", Marker="o", SizeData=350);
    scatter(XShu(1,t)/AU,XShu(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="black", Marker="^", SizeData=350);
    scatter(reshape(XSats(1,t,:), [Nsats,1])/AU,reshape(XSats(2,t,:), [Nsats,1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
    plot(reshape(XOutline(1,t,:)/AU, [outline_res, 1]), reshape(XOutline(2,t,:)/AU, [outline_res, 1]), LineStyle="--", LineWidth=2, Color='black')
    
    % Status indicator
    annotation('textbox', [0.8, 0.8, 0.175, 0.095], ...
                'String', stage(t), ...
                'FontSize', 20, ...
                'FontWeight', 'bold', ...
                'Color', 'white', ...
                'BackgroundColor', background(t), ...
                'EdgeColor', 'black', ...
                'HorizontalAlignment', 'center');

    % Add axis labels and title
    xlabel('x [AU]', 'FontWeight', 'bold');
    ylabel('y [AU]', 'FontWeight', 'bold');
    
    if t < end_of_deployment_idx
        title("Time since start of deployment procedure: " + num2str((etR(t)-etR_transfer(1))/(3600*24*365.25), '%.2f') + " Earth years")
    else

        title("Total deployment time: " + num2str((etR(end_of_deployment_idx)-etR_transfer(1))/(3600*24*365.25), '%.2f') + " Earth years")
    end

    
    hold off
    % pause(1e-2);

    % Export high-res frame
    filename = sprintf('temp_frames/frame_%04d.png', t);
    exportgraphics(fig, filename, Resolution=300, BackgroundColor='white');
    m = "Frames exported: " + num2str(t) + "/" + num2str(length(etR)) % Progress tracking
end