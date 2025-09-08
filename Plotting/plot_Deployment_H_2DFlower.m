%% Hohmann Insertion Plotting for 2D Flower Constellations

% Parameters

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;

outline_res = 500;

interval = 3600*24*1.0; % 1 days
sim_end = 3600*24*365.25*5; % 5 years
etR = 0:interval:sim_end; % Row vector of times between start and end date

KEa = [AU; 0; 0; 0; 0; 0];
XEa = propagateFromKeplerians(KEa, muSu, etR);

% Define Flower parameters to be represented

Nsats = 100;
M = 1;
W = 1; % For now, this should be negative... still working out some bugs on that end
aFLO = 1.4*AU;
eFLO = 0.3;
% aFLO = AU*abs(W/F)^(2/3);
% eFLO = abs(AU/aFLO - 1);

% Representation indices

shuttle_launch_idx = floor(length(etR)/16);

%% Generate shuttle trajectory

XShuttle_E = XEa(:,1:shuttle_launch_idx);

% Compute shuttle deployment orbit semi-major axis
aFD = aFLO*abs(M/W)^(2/3);

% Shuttle transfer orbit (circ 2 circ here)

at = 0.5*(AU+aFD);
f0 = 0;
et = (aFD-AU)/(aFD+AU);
[~,~,~,~,~,wt] = Cartesian2Keplerian(XEa(:,shuttle_launch_idx), muSu);
    
K_st_i = [at, et, 0, 0, wt, f0];
t_transfer = pi*sqrt(at^3 / muSu);
shuttle_arrival_idx = shuttle_launch_idx+1+ floor(t_transfer/interval);
etR_transfer = etR(shuttle_launch_idx+1:shuttle_arrival_idx);
    
XShuttle_transfer = propagateFromKeplerians(K_st_i, muSu, etR_transfer);


% Shuttle deployment orbit

etR_deploy = etR(shuttle_arrival_idx+1:end) * (-(W*M)/abs(W*M));
K_sd_i = [aFD; 0; 0; 0; 0; wt+pi];
XShuttle_deploy = propagateFromKeplerians(K_sd_i, muSu, etR_deploy);

XShu = cat(2, XShuttle_E, XShuttle_transfer, XShuttle_deploy);

%% Generate Flower Constellation Satellite Trajectories

% XSats = NaN([6, length(etR), Nsats]); % Sats default to not appearing
XSats = repmat(XShu, [1,1,Nsats]); % Sats default to being (visually) in the shuttle

TFLO = 2*pi*sqrt(aFLO^3/muSu);
deploy_interval = TFLO/Nsats * abs(M); % to test
deploy_interval_idxs = deploy_interval/interval;

rD = aFD;
rF = aFLO*(1 + eFLO);
at_s = (rD + rF)/2;
et_s = abs((rD-rF)/(rD+rF));

t_transfer_s = pi*sqrt(at_s^3 / muSu);
t_transfer_s_idxs = floor(t_transfer_s/interval);


for n = 1:Nsats
    % Transfer
    transfer_start_idx = shuttle_arrival_idx + 1 + floor((n-1)*deploy_interval_idxs);
    transfer_end_idx = transfer_start_idx + t_transfer_s_idxs;
    etR_transfer_s = etR(transfer_start_idx:transfer_end_idx);

    [~,~,~,~,~,wt_Shu] = Cartesian2Keplerian(XShu(:,transfer_start_idx), muSu);

    
    if aFD <= aFLO*(1+eFLO)
        wt_s = wt_Shu;
        f0_s = 0;
    else
        wt_s = wt_Shu+pi;
        f0_s = pi;
    end

    Ksn = [at_s, et_s, 0, 0, wt_s, f0_s];
    Kflo = [aFLO, eFLO, 0, 0 , wt_Shu, pi];
    XSats(:,transfer_start_idx:transfer_end_idx,n) = propagateFromKeplerians(Ksn, muSu, etR_transfer_s);
    XSats(:,transfer_end_idx:end,n) = propagateFromKeplerians(Kflo, muSu, etR(transfer_end_idx:end));
    % Final 2D Flower Orbit
    % if (W*F)/abs(W*F) < 0
    %     Kflo = [aFLO, eFLO, 0, 0 , wt_s, pi];
    % else
    %     Kflo = [aFLO, eFLO, 0, 0 , wt_s+pi, 0];
    % end
end

end_of_deployment_idx = transfer_end_idx;

%% Generate Flower Outline for Plotting (TODO)

% dtheta = pi*((1 + eFLO/2)^(3/2) - 1); % see paper
% KO = [aFLO, 0, 0, 0, 0, (K_sd_i(6)-dtheta)];
% XO = propagateFromKeplerians(KO, muSu, [0, -etR(shuttle_arrival_idx)]);
% XO = XO(:,2);
% XOutline = generateMOGNearPlanet(XO, muSu, etR, eFLO, 1, outline_res, false);


%% Plotting

fig_num = 3242342;
fig = figure(fig_num);
set(fig, 'Visible', 'off', 'Units', 'pixels', 'Position', [200 200 700 600], 'Resize', 'off');
% to avoid GUI rendering issues, turn visibility off

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

limval = ceil(max(max(max(XSats)))/AU); % AU
for t = 1:length(etR)
    clf(fig_num)

    scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
    
    hold on
    viscircles([0,0], aFD/AU, Color="black", LineWidth=1);
    viscircles([0,0], 1, Color="blue", LineWidth=1);

    scatter(XEa(1,t)/AU,XEa(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="blue", Marker="o", SizeData=350);
    scatter(XShu(1,t)/AU,XShu(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="black", Marker="^", SizeData=350);
    scatter(reshape(XSats(1,t,:), [Nsats,1])/AU,reshape(XSats(2,t,:), [Nsats,1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
    % plot(reshape(XOutline(1,t,:)/AU, [outline_res, 1]), reshape(XOutline(2,t,:)/AU, [outline_res, 1]), LineStyle="--", LineWidth=2, Color='black')
    
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