% 2D Flower Constellation Plotting/Debugging

% Parameters

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
day = 3600*24;
etR = 0:2*day:1500*day;
% etR = 1;

eMa = 0.0934;
aMa = 1.524*AU;

M = 2;
W = -3; % Could be better to flip the sign convention here out of convenience... todo?
E = 0;
% Nsats = ceil(100*max(abs(F),abs(W)));
Nsats = 500;
aF = 1.25*AU;
eF = 0.2;
% aF = AU*abs(W/M)^(2/3);
% eF = abs(AU/aF - 1);

% Ki = [aF, eF, 0, 0, -pi/2, 0];
% XSats = NSATSpropagateFromKeplerians_2DFlower(Ki,Nsats,W,F,etR,muSu);
% if (round(F) == F) && (round(W) == W)
%     XSats(:,:,Nsats+1) = XSats(:,:,1); % Close loop for plotting and analysis if this is a full 2DFlower
% end

% emin = eMa; % Mars
emin = 0.1; 
% emax = abs(AU/aF - 1);
emax = 0.6;
% Ki = [aF, (emin+emax)/2, 0, 0, pi/2, 0];
Ki = [aF, eF, 0, 0, pi/2, 0];
XSats = NSATSpropagateFromKeplerians_2DFlower_Es(Ki,Nsats,W,M,E,etR,muSu,emin,emax,"triangle");

% Plotting

limval = ceil(2*max(max(max(abs(XSats))))/AU + 0.01)/2; % AU
% Mars
KMa = [aMa; eMa; 0; 0; pi/2; 0];
XMa = propagateFromKeplerians(KMa, muSu, etR);

% Earth
KEa = [AU; 0; 0; 0; 0; pi/2];
XEa = propagateFromKeplerians(KEa, muSu, etR);

% Static Plot of the shape at t=1

figure(123)
scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
hold on
% plot(reshape(XSats(1,1,:), [size(XSats, 3),1])/AU, reshape(XSats(2,1,:), [size(XSats, 3),1])/AU, Color=[0.6 0 1]);
scatter(reshape(XSats(1,1,:), [size(XSats, 3),1])/AU,reshape(XSats(2,1,:), [size(XSats, 3),1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
% viscircles([0,0], 1, Color="blue", LineWidth=1);
% viscircles([0,0], 1.52, Color="red", LineWidth=1);
% viscircles([0,0], 0.72, Color="green", LineWidth=1);
xlabel("x [AU]");
ylabel("y [AU]");
hold off


%% Running Plot

fig_num = 1234;
fig = figure(fig_num);
set(fig, 'Visible', 'on', 'Units', 'pixels', 'Position', [200 200 700 600], 'Resize', 'off');
for t = 1:length(etR)
    clf(fig_num)
    scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
    
    hold on
    scatter(reshape(XSats(1,t,:), [size(XSats, 3),1])/AU,reshape(XSats(2,t,:), [size(XSats, 3),1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
    
    % Mars
    % plot(XMa(1,:)/AU,XMa(2,:)/AU,Color="red",LineWidth=1);
    % scatter(XMa(1,t)/AU,XMa(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="red", Marker="o", SizeData=350);
    
    % Earth
    viscircles([0,0], 1, Color="blue", LineWidth=1);
    scatter(XEa(1,t)/AU,XEa(2,t)/AU,MarkerEdgeColor="black", MarkerFaceColor="blue", Marker="o", SizeData=350);

    xlabel("x [AU]");
    ylabel("y [AU]");
    hold off
    pause(1e-2)

    % Export high-res frame
    % filename = sprintf('temp_frames/frame_%04d.png', t);
    % exportgraphics(fig, filename, Resolution=300, BackgroundColor='white');
    % m = "Frames exported: " + num2str(t) + "/" + num2str(length(etR)) % Progress tracking
end

