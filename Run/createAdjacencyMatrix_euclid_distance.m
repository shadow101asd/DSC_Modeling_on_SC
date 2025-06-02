function A = createAdjacencyMatrix_euclid_distance(Xs, cutoff)
%CREATEADJACENCYMATRIX_EUCLID_DISTANCE
%   Computes symmetric, time-varying adjacency matrix A of pairwise
%   Euclidean distances between nodes.
%
%   Inputs:
%       Xs    - [3, nT, N] positions of N nodes over nT timesteps
%       cutoff - optional [1 x nT] vector; distances above cutoff set to 0
%
%   Output:
%       A - [N, N, nT] adjacency matrix of distances

    [~, nT, N] = size(Xs);

    % Rearrange for vectorized distance calculation
    X = permute(Xs, [2, 3, 1]);         % [nT, N, 3]

    Xa = reshape(X, [nT, N, 1, 3]);     % [nT, N, 1, 3]
    Xb = reshape(X, [nT, 1, N, 3]);     % [nT, 1, N, 3]
    diffs = Xa - Xb;                    % [nT, N, N, 3]

    % Euclidean distances
    A = sqrt(sum(diffs.^2, 4));         % [nT, N, N]
    A = permute(A, [2, 3, 1]);          % [N, N, nT]

    % Optional cutoff masking
    if nargin == 2 && ~isempty(cutoff)
        cutoff_mat = reshape(cutoff, [1, 1, nT]);  % [1, 1, nT] for broadcasting
        A(A > cutoff_mat) = 0;
    end
end
