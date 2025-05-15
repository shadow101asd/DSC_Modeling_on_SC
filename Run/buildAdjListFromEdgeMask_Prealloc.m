function adjList = buildAdjListFromEdgeMask_Prealloc(endNodes, weights, mask, N)
% Fully preallocated adjacency list builder for masked graph

    E = endNodes(mask, :);
    W = weights(mask);

    if ~isnumeric(E)
        error('EndNodes must be numeric after masking');
    end

    adjList = cell(N, 1);

    % Precompute degrees
    degrees = accumarray([E(:,1); E(:,2)], 1, [N,1]);
    maxDeg = max(degrees);
    Aflat = zeros(N, maxDeg);
    counts = zeros(N, 1);

    for k = 1:size(E,1)
        u = E(k,1);
        v = E(k,2);
        if isinf(W(k)), continue; end  % Ignore Inf edges
        counts(u) = counts(u) + 1; Aflat(u, counts(u)) = v;
        counts(v) = counts(v) + 1; Aflat(v, counts(v)) = u;
    end

    for i = 1:N
        adjList{i} = Aflat(i, 1:counts(i));
    end
end
