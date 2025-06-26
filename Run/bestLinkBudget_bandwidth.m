function metric = bestLinkBudget_bandwidth(X1,X2,XSats,R)
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,bandwidths] = networkAnalysis_bandwidth(X1,X2,"Earth","Mars",XSats,R,0.9); % Pretty sure alpha needs to be less than 1 since our bandwidths are negative
    metric = -mean(bandwidths); % negative since ga minimizes the objective function
end