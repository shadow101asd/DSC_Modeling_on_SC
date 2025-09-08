function [] = plot_F2D_Static(fig_num, aFLO, eFLO, F, W, mu, plot_type, show_planet_circles)

% 2D Flower Constellation Plotting

% Parameters

AU = 149600000; % 1 AU in km
etR = 0;

resolution = ceil(1000*max(abs(F),abs(W)));

Ki = [aFLO, eFLO, 0, 0, pi/2, 0];
XSats = NSATSpropagateFromKeplerians_2DFlower(Ki,resolution,W,F,etR,mu);
if (round(F) == F) && (round(W) == W)
    XSats(:,:,resolution+1) = XSats(:,:,1); % Close loop for plotting and analysis if this is a full 2DFlower
end

% Plotting

limval = ceil(max(max(max(XSats)))/AU); % AU

% Static Plot of the Flower Constellation Shape

figure(fig_num)
clf(fig_num)
scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
hold on
if plot_type == "scatter"
    scatter(reshape(XSats(1,1,:), [size(XSats, 3),1])/AU,reshape(XSats(2,1,:), [size(XSats, 3),1])/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=200);
else
    plot(reshape(XSats(1,1,:), [size(XSats, 3),1])/AU, reshape(XSats(2,1,:), [size(XSats, 3),1])/AU, Color=[0.6 0 1]);
end

if show_planet_circles
    viscircles([0,0], 1, Color="blue", LineWidth=1);
    viscircles([0,0], 1.52, Color="red", LineWidth=1);
    viscircles([0,0], 0.72, Color="green", LineWidth=1);
end
hold off

xlabel("x [AU]");
ylabel("y [AU]");

end