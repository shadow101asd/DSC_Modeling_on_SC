function patterns = slots_1D_tile_nonperiodic_recursive(L,S)

    G = gcd(L,S);
    fs = factor(G);
    p_empty = zeros([1, L*S]);
    
    if fs == 1
        Lis = factor(L);

        if Lis == L
            Sis = factor(S);

            if Sis == S
                % We're in the base case, we only have the two trivial patterns
        
                p_clumps = p_empty;
                p_clumps(1:S) = 1;
                
                p_spaced = p_empty;
                p_spaced(1:L:end) = 1;
        
                patterns = unique(vertcat(p_clumps,p_spaced), "rows");
                return
            else
                si = max(Sis);
                % TODO
            end
        else
            li = max(Lis);
            % TODO
        end
    end

    % Otherwise, we make recursive calls using the largest gcd factor
    gi = max(fs);

    p_s = slots_1D_tile_nonperiodic_recursive(L,S/gi);
    p_R = repmat(p_s, [1, gi]);

    p_l = slots_1D_tile_nonperiodic_recursive(L/gi,S);
    p_P = cat(2, p_l, repmat(0*p_l, [1, gi-1]));

    patterns = unique(vertcat(p_R,p_P), "rows");
end