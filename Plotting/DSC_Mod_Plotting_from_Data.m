function f = DSC_Mod_Plotting_from_Data(figure_num, run_idx, Nf, Nl, view_in_3D, Planet_string, sat_config_name, t)
%DSC_Mod_Plotting_from_Data Summary of this function goes here
%   Detailed explanation goes here
    
    % Default behavior of the function
    if nargin <= 4 || isempty(view_in_3D)
        view_in_3D = false; 
    end
    if nargin <= 5 || isempty(Planet_string)
        Planet_string = 'Mars'; 
    end
    if nargin <= 6 || isempty(sat_config_name)
        sat_config_name = "RF_all8m";
    end

    % Paramaters
    AU = 1.496e8;

    % Load data from Data File
    data_file = "../Data/run"+run_idx+"/Nf"+int2str(Nf)+"Nl"+int2str(Nl)+".mat";
    load(data_file, "XSats_i", "bandwidths", "NSats", "X_opt")
    nsats = sum(NSats);
    
    % Load data from Input Files

    input_file = "../Inputs/run"+run_idx+".mat";
    load(input_file, "muSu", "XEa", "etR")
    PX_string = "X" + Planet_string(1:2);
    XP2 = load("../Inputs/run"+run_idx+".mat", PX_string).(PX_string);
    
    % Load and process satellite and terminal specs
    load("../Inputs/SE/"+sat_config_name+".mat", "Sat_specs", "Earth_specs", "Comms_specs");
    P2_specs = load("../Inputs/SE/"+sat_config_name+".mat", Planet_string + "_specs").(Planet_string + "_specs");
    P2_specs.pstring = Planet_string;

    % Process / recreate everything needed for Plotting

    if nargin >= 8 && ~isempty(t) % If specific frame is requested
        etR = [etR(1), etR(t)];
        XEa = [XEa(:,1), XEa(:,t)];
        XP2 = [XP2(:,1), XP2(:,t)];
    end

    % Recreate XSats from XSats_i
    XSats = propagateFromCartesians(XSats_i,muSu,etR);
    
    D_antennas = 8*ones([1,Nf]); % Assume all 8m for now
    R_1km = buildAdjacency_LinkBudgetMatrix(NSats, Comms_specs, Sat_specs, Earth_specs, P2_specs, D_antennas); % Units: bps*km^2
    
    alpha = 0.9;
    [graphslist,~,paths,Bandwidths] = networkAnalysis_bandwidth(XEa,XP2,"Earth",Planet_string,XSats,R_1km,alpha);

    
    Xs = collateXs(XEa,XP2,XSats);
    
    % Plotting
    nT = length(etR);
    color1 = "blue";
    color2 = "red";
    lim = ceil(2*max(max(max(abs(Xs))))/AU)/2;
    year = 365.25*3600*24; % 1 year in seconds

    for t=1:nT
        % Graph plot

        f = figure(figure_num);
        clf(figure_num)

        if view_in_3D
            set(gca, "XLim", [-lim, lim], "YLim", [-lim, lim]);
            view(45,45)
        else
            set(gca, "XLim", [-lim, lim], "YLim", [-lim, lim], "DataAspectRatio", [1, 1, 1]);
        end

        hold on

        % Prune graph to only keep the path edges
        path = paths{t};
        G = graphslist{t};
        % Consecutive pairs along the path
        s = path(1:end-1);
        ts = path(2:end);
        
        % Find the corresponding edge indices in G
        eidx = findedge(G, s, ts);      % returns 0 where an edge is missing
        eidx = eidx(eidx > 0);         % keep only existing edges
        
        % Keep only those edges (drop the rest)
        H = rmedge(G, setdiff(1:numedges(G), eidx));

        if view_in_3D
            % Plot graphs
            % Highlight Sun + Planets
            scatter3(0,0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black', LineWidth=0.75);
            scatter3(XEa(1,t)/AU,XEa(2,t)/AU,XEa(3,t)/AU, 2e2, 'o', color1,'filled','MarkerEdgeColor', 'black');
            scatter3(XP2(1,t)/AU,XP2(2,t)/AU,XP2(3,t)/AU, 2e2, 'o', color2,'filled','MarkerEdgeColor', 'black');
            
            % Highlight paths
            scatter3(reshape(Xs(1,t,:)/AU,nsats+2,1,1),reshape(Xs(2,t,:)/AU,nsats+2,1,1),reshape(Xs(3,t,:)/AU,nsats+2,1,1), ...
                MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=50);
            p1 = plot(H, 'XData', reshape(Xs(1,t,:)/AU,nsats+2,1,1), ...
                'YData', reshape(Xs(2,t,:)/AU,nsats+2,1,1), 'ZData', reshape(Xs(3,t,:)/AU,nsats+2,1,1), ...
                "EdgeColor", [0.3010 0.7450 0.9330]);
        else
            % Plot graphs
            % Highlight Sun + Planets
            scatter(0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black', LineWidth=0.75);
            scatter(XEa(1,t)/AU,XEa(2,t)/AU, 2e2, 'o', color1,'filled','MarkerEdgeColor', 'black');
            scatter(XP2(1,t)/AU,XP2(2,t)/AU, 2e2, 'o', color2,'filled','MarkerEdgeColor', 'black');
            
            % For legend entries...
            scatter(reshape(Xs(1,t,:)/AU,nsats+2,1,1),reshape(Xs(2,t,:)/AU,nsats+2,1,1), ...
                MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=75);
            plot(nan, nan, Color='red', LineWidth=6);
            
            % Actual graph plot
            p1 = plot(H, 'XData', reshape(Xs(1,t,:)/AU,nsats+2,1,1), ...
                'YData', reshape(Xs(2,t,:)/AU,nsats+2,1,1), "EdgeColor", 'red', ...
                NodeColor = [0.7, 0.7, 0.7], EdgeAlpha = 1.0, LineWidth=14);

            % Just so that they appear on top
            scatter(reshape(Xs(1,t,:)/AU,nsats+2,1,1),reshape(Xs(2,t,:)/AU,nsats+2,1,1), ...
                MarkerEdgeColor="black", MarkerFaceColor=[0.7,0.7,0.7], Marker="square", SizeData=75);
            scatter(XEa(1,t)/AU,XEa(2,t)/AU, 2e2, 'o', color1,'filled','MarkerEdgeColor', 'black');
            scatter(XP2(1,t)/AU,XP2(2,t)/AU, 2e2, 'o', color2,'filled','MarkerEdgeColor', 'black');
            
        end

        % Highlight paths
        % highlight(p1, paths{t}, 'LineWidth', 3, 'EdgeColor', color2);
        hold off
        % Add axis labels and title
        xlabel('x [AU]');
        ylabel('y [AU]');


        % Bandwidth Plot
        figure(figure_num+1)
        clf(figure_num+1)
        
        plot((etR-etR(1))/year, Bandwidths, LineWidth=1)
        hold on
        yline(mean(Bandwidths), LineWidth=2)
        yline(min(Bandwidths), Color="red", LineWidth=2)
        scatter((etR(t)-etR(1))/year, Bandwidths(t), MarkerFaceColor="red", SizeData=100);
        hold off
        xlabel("Simulation time [years]")
        ylabel("Bandwidth [Mbps]")
        legend("Network Bandwidth", "Mean Bandwidth", "Min Bandwidth")
        set(gca, FontSize=14)
        ylim([0 50*ceil(max(Bandwidths)/50)])
        
        pause(0.1);
   
    end

   

end

