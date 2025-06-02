function max_payload_mass = computeMaxPayloadMass(payload_mass_limit, DVs, m_dry, m_prop_initial, Isp, g0)
    % Make sure all units are in kg and m/s!

    eK = exp(sum(DVs)/(g0*Isp));
    m_payload = (m_prop_initial + m_dry*(1-eK))/(eK - 1); % To doublecheck algebra
    
    if m_payload >= 0
        max_payload_mass = min(m_payload, payload_mass_limit);
    else
        max_payload_mass = NaN; % Infeasible DV even without a payload
    end
end