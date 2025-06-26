function A = createAdjacencyMatrix_euclid_distance_block_dense(Xs, cutoff, blockSize)
%CREATEADJACENCYMATRIX_EUCLID_DISTANCE_BLOCK_DENSE
%   Computes symmetric, time-varying adjacency matrix of Euclidean distances
%   in blocks, using dense output matrix. Should have an identical output
%   to createAdjacencyMatrix_euclid_distance. 
%
%   Inputs:
%       Xs        - [3, nT, N] positions of N nodes over nT timesteps
%       cutoff    - [1 x nT] vector; distances above cutoff set to 0
%       blockSize - scalar; number of nodes to process per block
%
%   Output:
%       A - [N, N, nT] dense symmetric adjacency matrix

    [~, nT, N] = size(Xs);
    A = zeros(N, N, nT);  % Preallocate full result

    if nargin < 3 || isempty(blockSize)
        blockSize = 100;  % Default block size
    end

    for t = 1:nT
        Xt = squeeze(Xs(:, t, :))';  % [N x 3]
        cutoff_t = inf;
        if nargin >= 2 && ~isempty(cutoff)
            cutoff_t = cutoff(t);
        end

        for i1 = 1:blockSize:N
            i2 = min(i1 + blockSize - 1, N);
            Xi = Xt(i1:i2, :);  % [Bi x 3]

            for j1 = i1:blockSize:N
                j2 = min(j1 + blockSize - 1, N);
                Xj = Xt(j1:j2, :);  % [Bj x 3]

                Bi = i2 - i1 + 1;
                Bj = j2 - j1 + 1;

                Xi_rep = reshape(Xi, [Bi, 1, 3]);
                Xj_rep = reshape(Xj, [1, Bj, 3]);

                diffs = Xi_rep - Xj_rep;                    % [Bi x Bj x 3]
                dists = sqrt(sum(diffs.^2, 3));             % [Bi x Bj]

                % Apply cutoff if needed
                if isfinite(cutoff_t)
                    dists(dists > cutoff_t) = 0;
                end

                % Fill full matrix (both halves for symmetry)
                A(i1:i2, j1:j2, t) = dists;
                if i1 ~= j1
                    A(j1:j2, i1:i2, t) = dists';  % symmetry
                end
            end
        end
    end
end