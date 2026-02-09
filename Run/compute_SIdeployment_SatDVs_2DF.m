function [Sat_DVs, DPs_opt, pattern_opt, maxDV] = compute_SIdeployment_SatDVs_2DF(Nw,NM,aF,eF,NL,Npl,frac,resolution,mu)
    
    muSu = 1.327124400419393e+11;
    AU = 149600000; % 1 AU in km

    if nargin <= 8 || isempty(mu)
        mu = 1.327124400419393e+11; % Sun's gravitational parameter by default
    end
    if nargin <= 7 || isempty(resolution)
        resolution = 100; % Default gridsearch resolution
    end
    if nargin <= 6 || isempty(frac)
        frac = 1;
    end

    % grid_resolution = 10;
    
    patterns = slots_1D_tile_nonperiodic_recursive(NL,Npl); % Possible DP patterns from combinatorics problem

    DV_func = @(DP) computeOptimalPhasingDV_2DF_extended(Nw,-NM,eF,DP)*sqrt(mu/aF * AU/muSu); % scale correctly according to aF
    
    % P = size(patterns, 1);
    % best_DVs = NaN([1, P]);
    % P0_opts = best_DVs;
    % for p = 1:P
    %     pattern  = patterns(p,:);
    % 
    %     % Evaluate on the grid to get a good first guess
    %     P0_fg = Inf;
    %     best_idx = 0;
    %     P0s_grid = linspace(0, 2*pi, grid_resolution);
    %     for i = 1:grid_resolution
    %         if wrapper_func(P0s_grid(i), pattern, DV_func, frac) < P0_fg
    %             best_idx = i;
    %         end
    %     end
    % 
    %     % Run local optimizer near the known first guess
    %     options = optimset('Display','iter');
    %     [P0_opt, fval] = fminunc(@(P0) wrapper_func(P0, pattern, DV_func, frac), P0s_grid(best_idx), options)
    %     best_DVs(p) = fval;
    %     P0_opts(p) = P0_opt;
    % end

    DVs = computePhasingDVs_2DF_sequentialM(Nw, NM, eF, resolution, 1)*sqrt(mu/aF * AU/muSu); % scale correctly according to aF
    [DVs_pats, DPs_pats] = get_DVs_for_pattern(DVs, patterns, frac);

    [~, idx] = min(DVs_pats(:));
    [p_opt, idx_P0_opt] = ind2sub(size(DVs_pats), idx);
    
    pattern_opt = patterns(p_opt,:);

    % Run local optimizer near the best grid point
    options = optimset('Display','off'); % for debugging
    [P0_opt, ~] = fminunc(@(P0) wrapper_func(P0, pattern_opt, DV_func, frac), DPs_pats(idx_P0_opt), options);

    [maxDV, Sat_DVs, DPs_opt] = wrapper_func(P0_opt, pattern_opt, DV_func, frac);
end

function [DV, DVs, DPs] = wrapper_func(P0, pattern, DV_func, frac)
    NSats = sum(pattern);
    
    offsets = linspace(0,frac*2*pi,length(pattern)+1);
    offsets_pattern = pattern .* offsets(1:end-1);
    DPs = P0 + [0,nonzeros(offsets_pattern)']; % Keep only non-zero DV assignments corresponding to sat insertions
    DPs = wrapTo2Pi(DPs); % clean things up
    
    DVs = zeros([1, NSats]);
    for i = 1:NSats
        DVs(i) = DV_func(DPs(i));
    end
    
    DV = max(DVs);
end