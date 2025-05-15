function adjList = buildAdjListFromEdgeMask_FAST(G, mask)
    % Always use G.Edges to avoid bad input
    endNodes = G.Edges.EndNodes(mask,:);
    weights = G.Edges.Weight(mask);

    if ~isnumeric(endNodes)
        endNodes = reshape(findnode(G, endNodes(:)), [], 2);
    end

    N = numnodes(G);
    adjList = cell(N,1);

    % Preallocate neighbor lists
    degrees = accumarray([endNodes(:,1); endNodes(:,2)], 1, [N,1]);
    maxDeg = max(degrees);
    Aflat = zeros(N, maxDeg);
    counts = zeros(N, 1);

    for k = 1:size(endNodes,1)
        u = endNodes(k,1);
        v = endNodes(k,2);
        if isinf(weights(k)), continue; end
        counts(u) = counts(u) + 1; Aflat(u, counts(u)) = v;
        counts(v) = counts(v) + 1; Aflat(v, counts(v)) = u;
    end

    for i = 1:N
        adjList{i} = Aflat(i, 1:counts(i));
    end
end