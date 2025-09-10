function [meanB, minB, maxB, bandwidths] = bestLinkBudget_bandwidth(X1,X2,XSats,R,maxDR_Mbps,P2_name)
    if nargin <= 5
        P2_name = "Mars";
    end

    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,bandwidths] = networkAnalysis_bandwidth(X1,X2,"Earth",P2_name,XSats,R,0.9); % Pretty sure alpha needs to be less than 1 since our bandwidths are negative
    bandwidths = -bandwidths; % negative since ga minimizes the objective function

    if nargin >= 5
        % Then there is a cap on the feasible datarate passed in
        bandwidths = min(bandwidths, maxDR_Mbps);
    end

    meanB = mean(bandwidths); 
end