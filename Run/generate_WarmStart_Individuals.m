function Initial_Pop = generate_WarmStart_Individuals(Nmax, run_idx, Nf_max, Nl_max, empty_feature_snippet, Nf_real)

    if nargin <= 5
        Nf_real = Nf_max; % Default behavior if not passed in
    end
    
    
    % Collect Existing Datapoints

    Nf_range = 1:Nf_max;
    Nl_range = 1:(Nl_max+floor(Nmax/(Nf_max*2)));
    
    Bs = NaN([length(Nl_range), length(Nf_range)]);
    for Nf = Nf_range
        for Nl = Nl_range
            filename = "../Data/run"+run_idx+"/Nf"+int2str(Nf)+"Nl"+int2str(Nl)+".mat";
            
            if isfile(filename)
                % Extract and save data of interest
                Bs(Nl, Nf) = load(filename, "B").B;
            end
        end
    end

    % Find up to Nmax Best Existing Elements

    [vals, linearIndices] = maxk(Bs(:), Nmax);
    [Nls_selected, Nfs_selected] = ind2sub(size(Bs), linearIndices(~isnan(vals)));
    
    % Build InitialPop from Selected Genomes

    N_sel = length(Nls_selected);
    Initial_Pop = repmat(empty_feature_snippet, [N_sel, Nf_real]);
    
    for i = 1:N_sel
        X_opt = load("../Data/run"+run_idx+"/Nf"+int2str(Nfs_selected(i))+"Nl"+int2str(Nls_selected(i))+".mat").X_opt;
        
        % If backpropagating, adjust Nls of a random feature to make a valid
        % warm start. (in both directions actually)
        
        while true
            Nl_comp = sum(X_opt(2:length(empty_feature_snippet):end));
            if Nl_comp > Nl_max
                rand_idx_Nl = 2 + length(empty_feature_snippet)*(randi(Nfs_selected(i), 1, "double")-1);
                X_opt(rand_idx_Nl) = X_opt(rand_idx_Nl) - 1;
            elseif Nl_comp < Nl_max
                rand_idx_Nl = 2 + length(empty_feature_snippet)*(randi(Nfs_selected(i), 1, "double")-1);
                X_opt(rand_idx_Nl) = X_opt(rand_idx_Nl) + 1;
            else
                break
            end
        end
       

        Initial_Pop(i, 1:length(X_opt)) = X_opt;
    end

end