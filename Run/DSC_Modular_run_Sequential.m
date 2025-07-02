function [] = DSC_Modular_run_Sequential(Nl_range, Nf_range, run_idx, Planet_string)
    if nargin >= 4
        for Nl = Nl_range
            for Nf = Nf_range
                DSC_Modular_run(Nl, Nf, run_idx, Planet_string)
            end
        end
    else
        for Nl = Nl_range
            for Nf = Nf_range
                DSC_Modular_run(Nl, Nf, run_idx)
            end
        end
    end
end