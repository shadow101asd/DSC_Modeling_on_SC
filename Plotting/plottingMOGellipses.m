% Plotting to explain the MOG construction

AU = 149600000; % km
muSu = 1.327124400419393e+11;

aMOG = AU;
eMOG = 0.3;
XP = Keplerian2Cartesian(aMOG, eMOG, 0, 0, -pi/2, 0, muSu);

TMOG = 2*pi*sqrt(aMOG^3/muSu);
resolution = 500;
etR = linspace(0,TMOG, resolution);

NSats = 7;
outline_res = 1000;

XSats = generateMOGNearPlanet(XP, muSu, etR, eMOG, 1, NSats, false);
XOutline = generateMOGNearPlanet(XP, muSu, etR, eMOG, 1, outline_res, false);

% Plotting

appearance_idx = resolution/5;

figure(123123)

for t = 1:resolution
    clf(123123)
    s = scatter(0,0, "pentagram", "yellow", "filled", MarkerEdgeColor="black", SizeData=250, Linewidth=1.0);
    % legend("Sun", "AutoUpdate", "off")
    set(gca,'fontsize', 14)
    hold on
    for i = 1:NSats
        CO = colororder;
        col = CO(i,:);
        
        if t > appearance_idx || i == 1
            plot(XSats(1,:,i)/AU, XSats(2,:,i)/AU, LineWidth=2, Color=col);
            scatter(XSats(1,t,i)/AU, XSats(2,t,i)/AU, "filled", MarkerFaceColor=col, MarkerEdgeColor="black", SizeData=150)
        end
    end
    if t > appearance_idx
        plot(reshape(XOutline(1,t,:)/AU, [outline_res, 1]), reshape(XOutline(2,t,:)/AU, [outline_res, 1]), LineStyle="--", LineWidth=2, Color='black')   
    end
    hold off
    lim = 1.5;
    xlim([-lim lim])
    ylim([-lim lim])
    pbaspect([1 1 1])
    xlabel("x [AU]")
    ylabel("y [AU]")
    % exportgraphics(gcf,'MOGEllipses_pres.gif','Append',true);
    pause(1e-2)
end