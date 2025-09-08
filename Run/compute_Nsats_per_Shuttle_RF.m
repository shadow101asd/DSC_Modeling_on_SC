function [Nsats_per_shuttle, D_antennas, Npl_final, N_layers, active_constraint, m_persat] = compute_Nsats_per_Shuttle_RF(Shuttle_specs, Sat_specs, Npl, mass_func)
   % Note: This function was designed for Nedges = 2 or 1. Would need
   % modification for any larger number of antennas (in the RF case that is
   % driven by antenna mass and size constraints).

    Npl = floor(Npl); % Guarantee integrality

    % Extract satellite main bus sizing information
    a = Sat_specs.mainBusDims(1);
    b = Sat_specs.mainBusDims(2);
    h_bus = Sat_specs.mainBusDims(3);
    d_circum = sqrt(a^2 + b^2);

    D_payload = Shuttle_specs.maxPayloadDiameter;
    
    % Verify that the satellite bus fits in the payload bay
    if d_circum > D_payload
        Nsats_per_shuttle = 0;
        D_antennas = 0;
        return
    end

    % Compute Relay Satellite Antenna Diameter
    if Npl == 1
        D_antennas = D_payload;
        Npl_final = 1;
    else
        D_antennas = D_payload * sin(pi/Npl)/(1+sin(pi/Npl));

        if d_circum > D_antennas
            D_antennas = d_circum; % Is conservative
            Npl_final = floor(pi/asin(d_circum/(D_payload - d_circum)));
        else
            Npl_final = Npl;
        end
    end

    % Compute Antenna Depth (ignoring mount + fp apparatus for now)
    
    h_antenna = D_antennas / 16 / Sat_specs.FoverD;

    % Compute Feasible Number of Layers
    
    N_layers = floor(Shuttle_specs.maxPayloadHeightatFullDiameter/(Sat_specs.Nedges*h_antenna + h_bus));

    Nsats_V = N_layers * Npl_final;
    

    % Mass Considerations

    % m_antenna = Sat_specs.rho_antenna * pi*D_antennas^2/2 * sqrt(1/4 + 16*(Sat_specs.FoverD)^2); % This was incorrect.
    % m_antenna = Sat_specs.rho_antenna * pi*D_antennas^2/2 * sqrt(1/4 + 1/(16^2*(Sat_specs.FoverD)^2)); % This is what the formula should have been. Replaced with higher-fidelity equation
    
    % Replace with higher-fidelity equation (which is correct, and actually
    % yields smaller antenna masses! Need to eventually rerun Mercury +
    % Jupiter runs...)
    m_antenna = Sat_specs.rho_antenna * 64*pi/3 * (Sat_specs.FoverD)^2 * D_antennas * ((1/4 + 1/(64*(Sat_specs.FoverD)^2))^(3/2) - 1/8);

    m_persat = Sat_specs.mainBusMass + Sat_specs.Nedges * m_antenna;

    m_persat_V = mass_func(Nsats_V);
    m_persat_V = m_persat_V(1,1); % some functions have multiple outputs

    if m_persat_V >= m_persat
        active_constraint = "volume";
        Nsats_per_shuttle = Nsats_V;
    else
        active_constraint = "mass";
        
        % Let's just use an estimate for now - this is a lower bound on the
        % # of Satellites since phasing insertion procedures usually have
        % exponential fuel costs w.r.t. Nsats as an input parameter.
        Nsats_per_shuttle = floor(Nsats_V*m_persat_V/m_persat);
    end    
end