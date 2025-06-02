function excess_fuel_mass = computeFuelMargin(payload_M, DVs, m_dry, m_prop_initial, Isp, g0)
    % Make sure inputs are in kg, m/s, m/s2
    m0 = m_dry + m_prop_initial + payload_M;
    mf = m0 * exp(-sum(DVs)/(g0*Isp));

    excess_fuel_mass = mf - (m_dry + payload_M);
end