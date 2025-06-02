function [m_dry_per_sat, exitflag] = getSatDryMasses_CircularLoop_USA(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO)
    % Constants
    AU = 1.496e8; % 1 AU in km
    g0 = 9.81; % m/s^2
    
    % Find desired Shuttle orbital radius to ensure timely deployment

    Ru = Rf/(1+1/maxO)^(2/3);

    % Find desired v_inf (hyperbolic excess velocity)

    [~, ~, DV1, DV2] = computeHohmannTransferCirc2Circ(AU, Ru, mu);
    vinf = abs(DV1);
    
    % Compute Earth escape velocity from LEO
    muE = 3.986004e5; % km3/s2
    RE = 6378; % km
    a = 200; % LEO orbit altitude (km)
    v_escape = sqrt(2*muE/(RE + a)); % Escape velocity at LEO (a = 200 km) in km/s

    % Find desired DV for Earth escape + appropriate vinf
    V = sqrt(vinf^2 + v_escape^2);

    v_LEO = sqrt(muE/(RE+a));
    
    DV_LEO2EE = abs(V-v_LEO); % km/s
    
    DV_EE2CIRC = abs(DV2);

    total_ShuttleDV = (DV_LEO2EE + DV_EE2CIRC)*1e3; % m/s
    totalShuttlePayload_mass = computeMaxPayloadMass(maxShuttlePayload, total_ShuttleDV, Shuttle_dryMass, Shuttle_wetMass-Shuttle_dryMass, Shuttle_Isp, g0);
    
    if isnan(totalShuttlePayload_mass)
        m_dry_per_sat = NaN;
        exitflag = 0;
        return
    end

    % Satellite Deployment
    % Under this deployment scheme, all satellites have identical orbital
    % insertion DVs.

    sat_m0 = totalShuttlePayload_mass/NSats;

    [sat_DV, ~, ~, ~] = computeHohmannTransferCirc2Circ(Ru, Rf, mu); % km/s
    sat_DV = sat_DV*1e3; % Convert to m/s

    m_dry_per_sat = sat_m0*exp(-sat_DV/(Sat_Isp*g0));
    exitflag = 1; % Success!
end