function [Sats_MOGMs, Shuttle_FE, Sat_FEs, leftoverShuttle_FM, exitflag] = getMOG_SatDryMasses(aMOG,eMOG,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload,frac)
    
    if nargin <= 9
        frac = 1; % Default assumption is deploying a full MOG if no frac argument is passed in
    end

    [Sats_MOGMs_1, Shuttle_FE_1, Sat_FEs_1, leftoverShuttle_FM_1, exitflag_1] = getSatDryMasses2BI(aMOG,eMOG,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload);
    [Sats_MOGMs_2, Shuttle_FE_2, Sat_FEs_2, leftoverShuttle_FM_2, exitflag_2] = getSatDryMassesSI(aMOG,eMOG,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload,frac);
    
    % Sats_MOGMs_2 = 0;
    % exitflag_2 = 0; % Not considering MOG phasing at the moment

    if exitflag_1 == 0 && exitflag_2 == 0
        exitflag = 0;
        return
    end

    if exitflag_2 == 0 || Sats_MOGMs_1(1) >= Sats_MOGMs_2(1)
        Sats_MOGMs = Sats_MOGMs_1;
        Shuttle_FE = Shuttle_FE_1;
        Sat_FEs = Sat_FEs_1;
        leftoverShuttle_FM = leftoverShuttle_FM_1;
        exitflag = exitflag_1;
    else
        Sats_MOGMs = Sats_MOGMs_2;
        Shuttle_FE = Shuttle_FE_2;
        Sat_FEs = Sat_FEs_2;
        leftoverShuttle_FM = leftoverShuttle_FM_2;
        exitflag = exitflag_2;
    end
end