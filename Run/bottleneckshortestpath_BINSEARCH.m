function [path, D, edgepath] = bottleneckshortestpath_BINSEARCH(G, nodename1, nodename2, margin)
%BOTTLENECKSHORTESTPATH_BINSEARCH
    weights = G.Edges.Weight;

    idx1 = findnode(G, nodename1);
    idx2 = findnode(G, nodename2);
    
    % Init
    keep_idx = 1:length(weights);

    % if nargin == 4
    %     % Then we have a warm start
    %     D_prev = bottleneckDistance(G,previous_path);
    %     keep_idx = find(weights <= D_prev);
    % else
    %     % Use the minspantree() as a heuristic
    %     T = minspantree(G);
    %     D_heur = max(T.Edges.Weight);
    %     keep_idx = find(weights <= D_heur);
    % end

    while true
        lambda = mean(weights(keep_idx));
        
        % Check termination condition
        if abs(max(weights(keep_idx)) - min(weights(keep_idx))) < 1e-4
            D_bottleneck = max(weights(keep_idx));
            break
        end

        edges2remove = find(weights > lambda);

        Gprime = rmedge(G, edges2remove);

        bins = conncomp(Gprime);
        if bins(idx1) == bins(idx2)
            mask = weights(keep_idx) <= lambda;
        else
            mask = weights(keep_idx) > lambda;
        end
        keep_idx = keep_idx(mask);
    end

    % Reconstruct Path from Dequiv
    DMax = margin*D_bottleneck;
    
    % edges2remove_2 = find(weights > DMax);
    % G.Edges.OriginalIndex = (1:G.numedges)';
    % G2 = rmedge(G, edges2remove_2);
    % [path,D,new_edgepath] = shortestpath(G2,nodename1,nodename2);
    % edgepath = G2.Edges.OriginalIndex(new_edgepath);

    edges2remove_2 = find(weights > DMax);
    G.Edges.Weight(edges2remove_2) = Inf;
    [path,D,edgepath] = shortestpath(G,nodename1,nodename2);
    G.Edges.Weight = weights;
end

