function [m_persat, maxTOP, Sat_FEs] = computeSatMForPhasingDeployment(aMOG, eMOG, mu, NSats, Sat_Isp, Total_Payload_M, frac)
    g0 = 9.81; % m/s2
    ve = Sat_Isp*g0; % m/s
    
    if nargin <= 6
        frac = 1; % Default assumption is deploying a full MOG
    end
    [DVs, ToPs] = getDVsForPhasingDeployment(aMOG, eMOG, mu, NSats, frac); % outputs in [km/s, s]

    maxTOP = max(ToPs);
    DVs = DVs * 1000; % convert to m/s

    S = sum(exp(DVs/ve));
    m_persat = Total_Payload_M/S;
    Sat_FEs = m_persat*(exp(DVs/ve)-1); % Fuel mass expended by each satellite (kg)
end