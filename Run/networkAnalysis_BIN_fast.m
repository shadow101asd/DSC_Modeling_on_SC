function [graphslist, numedges, paths, Dpaths, Dequivs] = networkAnalysis_BIN_fast(X1, X2, name1, name2, XSats)
% High-efficiency network analysis with bottleneck path search
    margin = 1.01;

    % Create adjacency matrices
    Xs = collateXs(X1, X2, XSats);
    D12 = distanceBetweenXs(X1, X2);
    [~, nT, numsats] = size(XSats);
    A = createAdjacencyMatrix_euclid_distance(Xs(1:3,:,:), D12 * margin);

    % Create node names
    names = strings(numsats+2,1);
    for j = 1:numsats
        names(j) = "Sat" + j;
    end
    names(numsats+1) = name1;
    names(numsats+2) = name2;

    % Init outputs
    graphslist = cell(nT,1);
    numedges = zeros(nT,1);
    paths = cell(nT,1);
    Dpaths = zeros(nT,1);
    Dequivs = zeros(nT,1);

    for t = 1:nT
        At = A(:,:,t);
        G = graph(At, names);

        if t == 1
            [path, Dpath, edgepath] = bottleneckshortestpath_BINSEARCH_fast(G, name1, name2, margin, At);
        else
            D_prev = bottleneckDistance(G, paths{t-1});
            weights = G.Edges.Weight;
            keep_idx = find(weights <= D_prev * margin);
            [path, Dpath, edgepath] = bottleneckshortestpath_BINSEARCH_fast(G, name1, name2, margin, At, keep_idx);
        end

        graphslist{t} = G;
        numedges(t) = G.numedges;
        if isempty(path)
            error("No path found at timestep %d", t);
        end
        paths{t} = path;
        Dpaths(t) = Dpath;
        Dequivs(t) = bottleneckDistance2(G, edgepath);
    end
end
