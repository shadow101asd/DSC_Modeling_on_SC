function [DVs, DVs_seq] = computePhasingDVs_2DF_sequentialM(W, M, eF, resolution, Nm)
    dPs = linspace(0,2*pi,resolution);
    DVs = dPs*0;
    for i = 1:resolution
        DVs(i) = computeOptimalPhasingDV_2DF_extended(W, M, eF, dPs(i));
    end
    
    DVs_seq = findMultiPathOverPeriodicDistribution(DVs, Nm, resolution);
end