function [m_dry_per_sat, exitflag] = getSatDryMasses_CircularLoop_PA(Rf,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload, maxO)
    % Constants
    AU = 1.496e8; % 1 AU in km
    g0 = 9.81; % m/s^2

    % Shuttle Considerations
    % Find desired v_inf (hyperbolic excess velocity)

    [~, ~, DV1, DV2] = computeHohmannTransferCirc2Circ(AU, Rf, mu);
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

    % Satellite Phasing Considerations
    % We size the satellite wet masses such that their delivered "dry"
    % masses are all equal, similar to the Spread Out Insertion methodology
    % for MOG insertion.

    % Get phase angle changes (to streamline)

    if mod(NSats, 2) == 0
        DPhis = linspace(pi/NSats,2*pi*(1-1/(2*NSats)),NSats);
    else
        DPhis = linspace(0,2*pi*(1-1/NSats),NSats);
    end

    Dfs = DPhis/(2*pi);
    Dfs(ceil(NSats/2)+1:end) = abs(Dfs(ceil(NSats/2)+1:end) - 1);
    Dfs = Dfs/maxO; % Can achieve a given offset after maxO orbits

    % Get DVs
    
    a_pps = Rf*(1+Dfs).^(2/3);
    DVs = 2*abs(sqrt(mu/Rf)*(sqrt(2 - Rf./a_pps) - 1));
    DVs = DVs * 1000; % convert to m/s

    % Get dry mass per satellite

    S = sum(exp(DVs/(Sat_Isp*g0)));
    m_dry_per_sat = totalShuttlePayload_mass./S;

end