function [Sats_DryMasses, mode, exitflag, mpersat_under, mpersat_over, mpersat_phasing] = getSatDryMasses_CircularLoop(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO, frac)

% Compute Shuttle Orbit Options

mpersat_under = getSatDryMasses_CircularLoop_USA(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO);
mpersat_over = getSatDryMasses_CircularLoop_OSA(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO);
mpersat_phasing = getSatDryMasses_CircularLoop_PA(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO, frac);

ms = [mpersat_under, mpersat_over, mpersat_phasing];
[maxM, i] = max(ms);

if maxM > 0
    exitflag = 1;
    Sats_DryMasses = ms(i) * ones([NSats, 1]);
    if i == 1
        mode = "under";
    elseif i == 2
        mode = "over";
    else
        mode = "phasing";
    end
else
    exitflag = 0;
    Sats_DryMasses = [];
    mode = [];
end

end