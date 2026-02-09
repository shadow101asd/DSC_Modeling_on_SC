function [fig, p2] = plot_DSC_Modular_Structure(XSats_i, run_idx, figure_num, P2_name, SatSize)

    if nargin <= 3
        P2_name = 'Mars'; % Default behavior
    end

    if nargin <= 4
        SatSize = 200;
    end

    % Load Planetary Ephemerides from Inputs
    input_file = "../Inputs/run" + run_idx + ".mat";
    load(input_file, "XEa");
    XP2 = load(input_file, "X"+P2_name(1:2)).("X"+P2_name(1:2));

    % Plotting prep
    AU = 149600000; % 1 AU in km

    XSats = reshape(XSats_i, 6, []);
    XEa = XEa(:,1);
    XP2 = XP2(:,1);

    limval = ceil(max(max(abs(XSats)))/AU); % AU

    switch P2_name
        case 'Mars'
            P2_color = "red";
        case 'Venus'
            P2_color = "green";
        otherwise
            P2_color = "black";
    end

    % Plot
    fig = figure(figure_num);
    
    scatter(0,0, 450, 'pentagram', 'yellow', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    hold on
    scatter(XEa(1)/AU,XEa(2)/AU,MarkerEdgeColor="black", MarkerFaceColor="blue", Marker="o", SizeData=350);
    p2 = scatter(XP2(1)/AU,XP2(2)/AU,MarkerEdgeColor="black", MarkerFaceColor=P2_color, Marker="o", SizeData=350);
    sats = scatter(XSats(1,:)/AU,XSats(2,:)/AU,MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=SatSize);
    hold off
    xlabel("x [AU]", FontWeight="bold")
    ylabel("y [AU]", FontWeight="bold")

    set(fig, "Units", "pixels", "Position", [100 1200 600 600])
    set(gca, "XLim", [-limval, limval], "YLim", [-limval, limval], "DataAspectRatio", [1, 1, 1], "FontSize", 14);
end