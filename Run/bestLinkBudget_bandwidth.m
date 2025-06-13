function metric = bestLinkBudget_bandwidth(X1,X2,XSats,R)
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,bandwidths] = networkAnalysis_bandwidth(X1,X2,"Earth","Mars",XSats,R,1.2);
    metric = -mean(bandwidths); % negative since ga minimizes the objective function
end