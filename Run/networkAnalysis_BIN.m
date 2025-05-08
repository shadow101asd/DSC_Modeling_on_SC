function [graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis_BIN(X1,X2,name1,name2,XSats)
    % Network Analysis
    margin = 1.01;

    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats);
    D12 = distanceBetweenXs(X1,X2);
    [~,nT,numsats] = size(XSats);
    % A = createAdjacencyMatrix(Xs,@distanceBetweenXs, max(0.1*D12, 2*min(D12)));
    A = createAdjacencyMatrix2(Xs,@distanceBetweenXs, D12);

    names(numsats+1) = name1;
    names(numsats+2) = name2;
    for j=1:numsats
        names(j) = strcat("Sat", num2str(j));
    end

    % Initialize fields

    graphslist{nT} = [];
    numedges = zeros(nT,1);
    paths{nT} = {};
    Dpaths = zeros(nT,1);
    Dequivs = zeros(nT,1);
    
    for t = 1:nT
        At = A(:,:,t);
        
        g = graph(At, names);

        % Compute bottleneck shortest path
        if t == 1
            [path, Dpath, edgepath] = bottleneckshortestpath_BINSEARCH(g, name1, name2, margin);
        else
            D_prev = bottleneckDistance(g,paths{t-1});
            es2r = find(g.Edges.Weight > D_prev*margin);
            g = rmedge(g, es2r);
            [path, Dpath, edgepath] = bottleneckshortestpath_BINSEARCH(g, name1, name2, margin);
        end

        graphslist{t} = g;
        numedges(t) = g.numedges;

        if isempty(path)
            error("Failed at finding a path...")
        end
        
        paths{t} = path;
        Dpaths(t) = Dpath;
        Dbottleneck = bottleneckDistance2(g,edgepath);
        % Dbottleneck = bottleneckDistance(g,path);
        Dequivs(t) = Dbottleneck;
    end
end

