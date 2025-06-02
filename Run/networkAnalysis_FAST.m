function [graphslist,numedges,paths,Dpaths,Dequivs] = networkAnalysis_FAST(X1,X2,name1,name2,XSats,alpha)
    % Network Analysis

    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats);
    D12_i = distanceBetweenXs(X1(:,1),X2(:,1));
    [~,nT,numsats] = size(XSats);
    A = createAdjacencyMatrix_euclid_distance(Xs(1:3,:,:));

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
    
    t = 1;
    while t <= nT
        At = A(:,:,t);
        if t == 1
            At(At>D12_i) = 0.0;
        else
            At(At>alpha*Dequivs(t-1)) = 0.0;
        end
        
        g = graph(At, names);
        graphslist{t} = g;
        numedges(t) = g.numedges;

        % Compute bottleneck shortest path
        [path, Dpath, ~, DBottleneck] = bottleneckshortestpath3(g, name1, name2, false);
        
        % Catch errors
        if isempty(path) 
            % Adjust alpha and try again
            alpha = alpha + 0.1;
            continue
        else % We're good to go
            paths{t} = path;
            Dpaths(t) = Dpath;
            Dequivs(t) = DBottleneck;
            t = t + 1; % Move on to next loop
        end
        
    end
end

% TODO: Implement following along previous path as a more rigorous upper
% bound?