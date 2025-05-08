function f = DSC_TVG_Plotting_from_Xs_2P(figure_num, XSats, X1, X2, X3, name1, name2, name3, color1, color2, color3, etR)
%DSC_TVG_Plotting Summary of this function goes here
%   Detailed explanation goes here
    AU = 1.496e8;

    X12s = collateXs(X1,X2,XSats);
    X13s = collateXs(X1,X3,XSats);
    [graphslist12,~,paths12,~,Dequivs12] = networkAnalysis(X1,X2,name1,name2,XSats);
    [graphslist13,~,paths13,~,Dequivs13] = networkAnalysis(X1,X3,name1,name3,XSats);
    [~,nT] = size(X1);
    [~,~,nsats] = size(XSats);
    
    % Plotting
    lim = 1.5; % AU
    f = figure(figure_num);
    for t=1:nT
        clf(figure_num)
        % viscircles([0,0], 2.2, Color="black");
        set(gca, "XLim", [-lim, lim], "YLim", [-lim lim], "DataAspectRatio", [1, 1, 1]);
        hold on

        % Plot paths
        g12 = graphFromPath(graphslist12{t}, paths12{t});
        g13 = graphFromPath(graphslist13{t}, paths13{t});
        p1 = plot(g12, 'XData', reshape(X12s(1,t,:)/AU,nsats+2,1,1), ...
            'YData', reshape(X12s(2,t,:)/AU,nsats+2,1,1), "EdgeColor", color2, ...
            "NodeColor", [0.3 0.60 1.0], "LineWidth", 7);
        p2 = plot(g13, 'XData', reshape(X13s(1,t,:)/AU,nsats+2,1,1), ...
            'YData', reshape(X13s(2,t,:)/AU,nsats+2,1,1), "EdgeColor", color3, ...
            "NodeColor", [0.3 0.60 1.0], "LineWidth", 7);

        % Highlight Sun + Planets
        scatter(0,0, 5e2, 'pentagram', 'yellow','filled', 'MarkerEdgeColor', 'black');
        scatter(X1(1,t)/AU,X1(2,t)/AU,2e2, 'o', color1,'filled','MarkerEdgeColor', 'black');
        scatter(X2(1,t)/AU,X2(2,t)/AU,2e2, 'o', color2,'filled','MarkerEdgeColor', 'black');
        scatter(X3(1,t)/AU,X3(2,t)/AU,2e2, 'o', color3,'filled','MarkerEdgeColor', 'black');

        hold off
        % Add axis labels and title
        xlabel('x [AU]');
        ylabel('y [AU]');
        title(sprintf(['ICONIC Architecture B, %d Satellites. \n ' ...
            '%s-%s: Resulting ED: %f AU \n %s-%s: Resulting ED: %f AU'], ...
            nsats, name1, name2, mean(Dequivs12)/AU, name1, name3, mean(Dequivs13)/AU));
        pause(0.1);
    end
    
    function g = graphFromPath(gold, path)
        A = false(height(gold.Nodes));
        g = graph(A, gold.Nodes);
        for n = 1:(length(path)-1)
            g = addedge(g, path(n), path(n+1));
        end
    end
end

