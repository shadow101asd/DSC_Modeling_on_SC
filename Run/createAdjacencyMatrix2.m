function A = createAdjacencyMatrix2(Xs, pairwiseFunc, cutoff)
    % Create a symmetric, time-varying adjacency matrix A of size [N,N,nT]
    % Xs is [dim, nT, N] — dim = spatial dimensions, nT = #timesteps, N = #nodes
    % pairwiseFunc(Xa, Xb) must return a [1, nT] vector of weights
    % cutoff is optional, [1 x nT]

    [~, nT, N] = size(Xs);
    A = zeros(N, N, nT);

    % Only compute upper triangle (a < b)
    for a = 1:N-1
        Xa = Xs(:, :, a);  % Pre-extract
        for b = a+1:N
            Xb = Xs(:, :, b);  % Pre-extract
            w = pairwiseFunc(Xa, Xb);  % Should return [1, nT]
            A(a, b, :) = w;
            A(b, a, :) = w;  % Symmetric
        end
    end

    % Optional cutoff masking
    if nargin == 3 && ~isempty(cutoff)
        for t = 1:nT
            At = A(:, :, t);
            At(At > cutoff(t)) = 0.0;
            A(:, :, t) = At;
        end
    end
end