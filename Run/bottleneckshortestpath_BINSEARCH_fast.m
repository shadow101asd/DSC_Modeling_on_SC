function [path, D, edgepath] = bottleneckshortestpath_BINSEARCH_fast(G, nodename1, nodename2, margin, At, keep_idx)
% FINAL FAST VERSION: weight masking + shortestpath
% MATLAB-optimized: no rmedge, no adjacency lists

    weights = G.Edges.Weight;
    idx1 = findnode(G, nodename1);
    idx2 = findnode(G, nodename2);

    if nargin <= 5
        keep_idx = 1:length(weights);
    end

    while true
        w = weights(keep_idx);
        lambda = mean(w);

        if abs(max(w) - min(w)) < 1e-4
            D_bottleneck = max(w);
            break
        end

        At2 = At;
        At2(At > lambda) = 0;  

        % Check connectivity
        if bfs_reachable_At(At2, idx1, idx2)
            keep_idx = keep_idx(w <= lambda);
        else
            keep_idx = keep_idx(w > lambda);
        end
    end

    % Final step: shortestpath with margin-masked edges
    finalMask = weights <= (margin * D_bottleneck);
    G.Edges.Weight(~finalMask) = Inf;

    [path, D, edgepath] = shortestpath(G, nodename1, nodename2);
end