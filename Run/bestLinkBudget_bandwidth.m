function metric = bestLinkBudget_bandwidth(X1,X2,XSats,R,maxDR_Mbps)
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,bandwidths] = networkAnalysis_bandwidth(X1,X2,"Earth","Mars",XSats,R,0.9); % Pretty sure alpha needs to be less than 1 since our bandwidths are negative

    if nargin >= 5
        % Then there is a cap on the feasible datarate passed in
        bandwidths = min(bandwidths, maxDR_Mbps);
    end

    metric = -mean(bandwidths); % negative since ga minimizes the objective function
end