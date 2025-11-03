function [DVs_for_patterns, desired_grid] = get_DVs_for_pattern(DVs, patterns, frac)
    
    nV = length(DVs);
    resolution = max(100, nV);
    interval = 2*pi*frac/size(patterns,2); % DP interval between adjacent pattern slots
    N = floor(2*pi/interval);
    div = ceil(resolution/N); % Should be at least 1 this way
    desired_grid = 0:interval/div:(2*pi); % grid of DPs
    
    DPs = linspace(0, 2*pi, nV);
    interpolated_DVs = interp1(DPs,DVs,desired_grid);

    % For each pattern, do the thing
    patterns_spaced = zeros([size(patterns,1), size(patterns,2)*div]);
    patterns_spaced(:,1:div:end) = patterns;
    DVs_for_patterns = zeros(size(patterns,1),length(interpolated_DVs));

    for p = 1:size(patterns,1)
        for i = 1:length(interpolated_DVs)
            idxs = mod(find(patterns_spaced(p,:))+i-1,length(interpolated_DVs));
            idxs(idxs==0) = length(interpolated_DVs); % clean up mod output
            DVs_for_patterns(p,i) = max(interpolated_DVs(idxs));
        end
    end

end