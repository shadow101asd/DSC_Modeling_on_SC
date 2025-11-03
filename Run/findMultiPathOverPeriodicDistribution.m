function seq = findMultiPathOverPeriodicDistribution(raw, N, period)
% findMultiPathOverPeriodicDistribution
% seq(r) = min over i1,...,iN of raw(i1)+...+raw(iN)
%          such that (i1+...+iN) ≡ r (mod period)
%
% Notes:
% - 1-based indexing is used (MATLAB). Residues are taken in 1..period.
% - Reuse of indices is allowed (if you need indices to be distinct,
%   this requires a different DP).
%
% Time complexity: O(N * period^2)

    if nargin < 3 || isempty(period)
        period = numel(raw);
    end
    P = period;
    L = numel(raw);
    raw = raw(:).';  % row vector

    % Compress raw into minimal cost per residue class modulo P.
    % base(r) = min raw(i) over all i with i ≡ r (mod P)
    base = inf(1, P);
    for i = 1:L
        r = mod(i-1, P) + 1;     % residue in 1..P
        if raw(i) < base(r)
            base(r) = raw(i);
        end
    end

    % If N <= 1: return the original series when P == L (natural case),
    % otherwise return the residue-compressed base.
    if N <= 1
        if P == L
            seq = raw(1:P);
        else
            seq = base;
        end
        return
    end

    % DP over the min-plus semiring with cyclic indices:
    % dp_n = dp_{n-1} (⊗) base, where (⊗) is min-plus cyclic convolution.
    dp = base;
    for step = 2:N
        newdp = inf(1, P);
        % newdp(r) = min_k [ dp(k) + base( (r - k) mod P ) ]
        for r = 1:P
            best = inf;
            for k = 1:P
                idx = mod(r - k, P) + 1;  % cyclic shift
                val = dp(k) + base(idx);
                if val < best
                    best = val;
                end
            end
            newdp(r) = best;
        end
        dp = newdp;
    end

    seq = dp;
end
