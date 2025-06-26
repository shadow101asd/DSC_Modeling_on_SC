function [] = DSC_Modular_run(Nl, Nf, run_idx)
    gaoptions = optimoptions("ga", Display="iter", FunctionTolerance=1e-5, ...
        MaxStallGenerations=20, PopulationSize=300, UseParallel=false);

    sat_config_name = "RF_all8m";
    shuttle_config_name = "SS_BIII";

    % Run
    DSC_Modular_RF_2P(Nl, Nf, run_idx, gaoptions, sat_config_name, shuttle_config_name, false)
end