function [] = DSC_Modular_run(Nl, Nf, run_idx, Planet_string)
    gaoptions = optimoptions("ga", Display="iter", FunctionTolerance=1e-5, ...
        MaxStallGenerations=3, PopulationSize=20, UseParallel=true);
    

    sat_config_name = "RF_all8m";
    shuttle_config_name = "SS_BIII";
    plotting_bool = true;
    warm_start_bool = true;

    % Run
    if nargin >= 4
        DSC_Modular_RF_2P(Nl, Nf, run_idx, gaoptions, sat_config_name, shuttle_config_name, plotting_bool, warm_start_bool, Planet_string)
    else
       DSC_Modular_RF_2P(Nl, Nf, run_idx, gaoptions, sat_config_name, shuttle_config_name, plotting_bool, warm_start_bool)
    end
end