function [path, D, edgepath] = mst_path(T, s, t)
    if ischar(s) || isstring(s), s = findnode(T, s); end
    if ischar(t) || isstring(t), t = findnode(T, t); end

    adjList = buildAdjListFromAdjacency(T);

    N = numnodes(T);
    visited = false(N, 1);
    parent = zeros(N, 1);
    queue = zeros(N, 1);
    front = 1; rear = 1;
    queue(rear) = s;
    visited(s) = true;

    while front <= rear
        u = queue(front); front = front + 1;
        nbrs = adjList{u};  % FAST

        for i = 1:length(nbrs)
            v = nbrs(i);
            if ~visited(v)
                visited(v) = true;
                parent(v) = u;
                rear = rear + 1;
                queue(rear) = v;
                if v == t, break; end
            end
        end
        if visited(t), break; end
    end

    if ~visited(t)
        path = [];
        D = Inf;
        edgepath = [];
        return;
    end

    % Reconstruct path
    path = zeros(1, N); k = 1; u = t;
    while u ~= s
        path(k) = u;
        u = parent(u);
        k = k + 1;
    end
    path(k) = s;
    path = path(k:-1:1);

    edgepath = findedge(T, path(1:end-1), path(2:end));
    D = sum(T.Edges.Weight(edgepath));

end

function adjList = buildAdjList(G)
    N = numnodes(G);
    E = G.Edges.EndNodes;

    % Convert to numeric node indices if needed
    if ~isnumeric(E)
        nodeIDs = findnode(G, E);  % One call
        E = reshape(nodeIDs, [], 2);  % Back to M×2
    end

    M = size(E, 1);  % number of edges

    % First pass: count neighbors
    deg = accumarray([E(:,1); E(:,2)], 1, [N, 1]);

    % Preallocate adjacency storage
    maxDeg = max(deg);
    A = zeros(N, maxDeg);
    counts = zeros(N, 1);

    % Second pass: fill
    for k = 1:M
        u = E(k,1); v = E(k,2);
        counts(u) = counts(u) + 1; A(u, counts(u)) = v;
        counts(v) = counts(v) + 1; A(v, counts(v)) = u;
    end

    % Convert to final cell array
    adjList = cell(N, 1);
    for i = 1:N
        adjList{i} = A(i, 1:counts(i));
    end
end

function adjList = buildAdjListFromAdjacency(G)
% Efficiently builds adjacency list from sparse adjacency matrix of G
% Works fastest with numeric node IDs

    A = adjacency(G);  % sparse matrix, N×N
    [i, j] = find(triu(A, 1));  % only upper triangle (i < j) for undirected edges
    N = numnodes(G);
    adjList = cell(N, 1);

    % Preallocate neighbor counts
    deg = accumarray([i; j], 1, [N, 1]);
    maxDeg = max(deg);
    Aflat = zeros(N, maxDeg);
    counts = zeros(N, 1);

    % Fill neighbors without duplicates
    for k = 1:length(i)
        u = i(k); v = j(k);

        counts(u) = counts(u) + 1;
        Aflat(u, counts(u)) = v;

        counts(v) = counts(v) + 1;
        Aflat(v, counts(v)) = u;
    end

    % Convert to cell array
    for n = 1:N
        adjList{n} = Aflat(n, 1:counts(n));
    end
end

