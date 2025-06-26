function metric = bestLinkBudget_bandwidth_NPs(X1,XPs,XSats,R)
    % TLDR: What's the metric that we're optimizing over?
    if size(XPs,3) == 1 % we return to the usual 2P scenario
        metric = bestLinkBudget_bandwidth(X1,XPs,XSats,R);
    else
        bandwidths = networkAnalysis_maxminB(X1,XPs,"Earth",["Mars","Venus"],XSats,R);
        metric = -mean(bandwidths); % negative since ga minimizes the objective function
    end
end