function patterns = slots_1D_tile_nonperiodic_recursive(L,S)

    p_empty = zeros([1, L*S]);
    
    Lis = factor(L);
    Sis = factor(S);

    if ~isscalar(Lis)
        li = max(Lis);
        p_l = slots_1D_tile_nonperiodic_recursive(L/li,S);
        p_l1 = cat(2, p_l, repmat(0*p_l, [1, li-1]));
        p_l2 = insert_zero_columns(p_l, 1, li-1);
        patterns_L = unique(vertcat(p_l1,p_l2), "rows");
    else
        patterns_L = [];
    end

    if ~isscalar(Sis)
        si = max(Sis);
        p_s = slots_1D_tile_nonperiodic_recursive(L,S/si);
        p_s1 = repmat(p_s, [1, si]);
        p_s2 = repelem(p_s, 1, si);
        patterns_S = unique(vertcat(p_s1,p_s2), "rows");
    else
        patterns_S = [];
    end

    % We always have the two trivial patterns

    p_clumps = p_empty;
    p_clumps(1:S) = 1;

    p_spaced = p_empty;
    p_spaced(1:L:end) = 1;

    patterns_base = unique(vertcat(p_clumps,p_spaced), "rows");

    % Concatenate the patterns from the recursive calls earlier
    patterns = unique(vertcat(patterns_S,patterns_L,patterns_base), "rows");

end

function B = insert_zero_columns(A, L, Y)
% B = insert_zero_columns(A, L, Y)
% Insert Y columns of zeros after every L columns of A,
% including after the last block if size(A,2) is a multiple of L.
%
% Example:
%   A = reshape(1:12,3,4);
%   B = insert_zero_columns(A, 2, 1);
%   % Adds one zero column after every 2 columns, including at the end if needed.

[nRows, nCols] = size(A);
nBlocks = ceil(nCols / L);
B = [];

for i = 1:nBlocks
    % Take the next L columns (or fewer for the last block)
    startIdx = (i-1)*L + 1;
    endIdx = min(i*L, nCols);
    Ablock = A(:, startIdx:endIdx);

    % Append block and zeros
    B = [B, Ablock, zeros(nRows, Y)];
end

end
