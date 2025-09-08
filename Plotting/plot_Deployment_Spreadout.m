%% Spread-out Insertion Plotting for AAS/AIAA MOG Presentation

% Parameters

AU = 149600000; % 1 AU in km

outline_res = 500;

interval = 3600*24*0.5; % 1 day
sim_end = 3600*24*365.25*5; % 5 years
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

% XSats = NaN([6, length(etR), Nsats]); % Sats default to not appearing
XSats = repmat(XShu, [1,1,Nsats]); % Sats default to being in the shuttle

% Compute Phase Changes

spacing = 2*pi/Nsats;
if rem(Nsats, 2) == 0 % If NSats is even
    dPs = linspace(spacing/2, 2*pi-spacing/2, Nsats);
else % If NSats is odd
    dPs = linspace(0, 2*pi - spacing, Nsats);
end

transfer_start_idxs = ones([Nsats 1]);
transfer_end_idxs = ones([Nsats 1]);
deployment_start_idx = shuttle_arrival_idx+1;

alpha = wt + pi/2; % Angle to rotate velocities by
Rz = [cos(alpha), -sin(alpha), 0;
      sin(alpha),  cos(alpha), 0;
           0    ,      0     , 1]; % Resulting rotation matrix

for n = 1:Nsats
    % Transfer
    [~, ~, ~, ~, details] = computeOptimalPhasingDVMOG_extended(aMOG, eMOG, dPs(n), muSu);
    transfer_start_idxs(n) = deployment_start_idx + round(details.t1/interval);
    transfer_end_idxs(n) = deployment_start_idx + round(details.t2/interval);
    V1 = Rz * (details.V1)'; % rotate to appropriate frame
    R1 = Rz * (details.R1)'; % rotate to appropriate frame

    % Xtransfer = [XShu(1:3, transfer_start_idxs(n)); V1];
    Xtransfer = [R1; V1];
    [a_ts, e_ts, i_ts, Om_ts, w_ts, f0_ts] = Cartesian2Keplerian(Xtransfer, muSu);
    Ktransfer = [a_ts, e_ts, i_ts, 0, w_ts, f0_ts];
    XSats(:,transfer_start_idxs(n):transfer_end_idxs(n),n) = propagateFromKeplerians(Ktransfer,muSu,etR(transfer_start_idxs(n):transfer_end_idxs(n)));
    
    % MOG Orbit: we know what the final orbit should be
    [~, ~, ~, ~, ~, f0_EOT] = Cartesian2Keplerian(XShu(:,transfer_end_idxs(n)), muSu);
    Kshu_at_EOT = [aMOG, eMOG, 0, 0, wt, f0_EOT];
    KMOG = shiftK_inMOG(Kshu_at_EOT, dPs(n), muSu);
    XSats(:,(transfer_end_idxs(n)):end,n) = propagateFromKeplerians(KMOG,muSu,etR((transfer_end_idxs(n)):end));
end

end_of_deployment_idx = max(transfer_end_idxs);

%% Generate MOG Outline for Plotting

XO = propagateFromKeplerians([aMOG, 0, 0, 0, 0, K_sd_i(5)+pi], muSu, [0, -etR(shuttle_arrival_idx)]);
XO = XO(:,2);
XOutline = generateMOGNearPlanet(XO, muSu, etR, eMOG, 1, outline_res, false);


%% Plotting

fig_num = 3242344;
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
    plotEllipse(K_sd_i, muSu, 1, AU, "black")
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
    pause(1e-2);

    % Export high-res frame
    filename = sprintf('temp_frames/frame_%04d.png', t);
    exportgraphics(fig, filename, Resolution=300, BackgroundColor='white');
    m = "Frames exported: " + num2str(t) + "/" + num2str(length(etR)) % Progress tracking
end


%% Assemble frames

v = VideoWriter('PlottingData/AAS Paper Figs and Data/MOG_Deployment_SI_pres2.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 100;
open(v);


for t = 1:length(etR)
    filename = sprintf('temp_frames/frame_%04d.png', t);
    img = imread(filename);
    img = im2uint8(img);
    writeVideo(v, img);
    m = "Frames assembled: " + num2str(t) + "/" + num2str(length(etR)) % Progress tracking
end


close(v);