function [Bandwidths, paths] = networkAnalysis_maxminB(X1,XPs,name1,namePs,XSats,R)
    % For now, assumes Nedges = 2 for all satellites, i.e. all paths are
    % exclusive.

    % Create Adjacency Matrices
    Xs = collateXs(X1,XPs,XSats);
     
    if isempty(XSats)
        numsats = 0;
        [~,nT,~] = size(X1); % more robust
    else
        [~,nT,numsats] = size(XSats);
    end
    
    A = createAdjacencyMatrix_euclid_distance_block_dense(Xs(1:3,:,:)); % block version used for avoiding memory issues for large Nsats
    B = -R ./ A.^2 * 1e-6; % Bandwidth adjacency matrix, Mbps. Negative signs are because bsp3 attempts to minimize edge weights.
    
    % Run integer linear program over each timestep
    
    x0 = [];
    Bandwidths = zeros([nT,1]);

    for t = 1:nT
        Bt = B();
        [Band,x0] = find_maxminB(Bt, x0, numsats);
        Bandwidths(t) = Band;
    end
end