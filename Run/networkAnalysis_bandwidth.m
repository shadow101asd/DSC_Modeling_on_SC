function [graphslist,numedges,paths,Bandwidths] = networkAnalysis_bandwidth(X1,X2,name1,name2,XSats,R,alpha)
    % Network Analysis - based on networkAnalysis_FAST, but with bandwidths
    % instead.

    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats);
    D12_i = distanceBetweenXs(X1(:,1),X2(:,1));
    B12_i = -R(end, end-1) / D12_i^2 * 1e-6; % Initial DTE bandwidth, Mbps. Negative signs are because bsp3 attempts to minimize edge weights.
    
    if isempty(XSats)
        numsats = 0;
        [~,nT,~] = size(X1); % more robust
    else
        [~,nT,numsats] = size(XSats);
    end
    
    A = createAdjacencyMatrix_euclid_distance_block_dense(Xs(1:3,:,:)); % block version used for avoiding memory issues for large Nsats
    B = -R ./ A.^2 * 1e-6; % Bandwidth adjacency matrix, Mbps. Negative signs are because bsp3 attempts to minimize edge weights.

    names(numsats+1) = name1;
    names(numsats+2) = name2;
    for j=1:numsats
        names(j) = strcat("Sat", num2str(j));
    end

    % Initialize fields

    graphslist{nT} = [];
    numedges = zeros(nT,1);
    paths{nT} = {};
    Bandwidths = zeros(nT,1);
    
    t = 1;
    while t <= nT
        Bt = B(:,:,t);
        Bt(1:size(Bt,1)+1:end) = 0; % Set self-loops to zero bandwidth
        Bt(isnan(Bt)) = 0; % make sure there are no NaNs left.
        Bt = min(Bt, Bt'); % Enforce symmetry (for now), taking largest bandwidth (most negative) option each time.

        if t == 1
            Bt(Bt>B12_i) = 0.0;
        else
            Bt(Bt>-alpha*Bandwidths(t-1)) = 0.0; % Needs to be negative since Bandwidths are positive and Bt are negative
        end
        
        g = graph(Bt, names, 'upper'); % Upper adds another redundancy to deal with weight symmetry issues
        graphslist{t} = g;
        numedges(t) = g.numedges;

        % Compute bottleneck shortest path
        [path, ~, ~, DBottleneck] = bottleneckshortestpath3(g, name1, name2, false);
        
        % Catch errors
        if isempty(path) 
            % Adjust alpha and try again. This works regardless of whether
            % alpha is larger or smaller than 1.
            alpha = alpha^2;
            continue
        else % We're good to go
            paths{t} = path;
            Bandwidths(t) = -DBottleneck; % Have to invert the negative sign again
            t = t + 1; % Move on to next loop
        end
        
    end
end

% TODO: Implement following along previous path as a more rigorous upper
% bound?