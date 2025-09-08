function [Sats_MOGMs, Shuttle_FE, Sat_FEs, leftoverShuttle_FM, exitflag] = getSatDryMassesSI(aMOG,eMOG,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload,frac)
%getSatDryMassesSI Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 9
    frac = 1; % Default assumption is deploying a full MOG
end

% Constants
AU = 1.496e8; % 1 AU in km
g0 = 9.81; % m/s^2

% Find desired v_inf
[~, ~, DV1, DV2] = computeHohmannTransferCirc2Ell(AU, aMOG, eMOG, mu); % outputs in km/s
vinf = abs(DV1);

muE = 3.986004e5; % km3/s2
RE = 6378; % km
a = 200; % LEO orbit altitude (km)
v_escape = sqrt(2*muE/(RE + a)); % Escape velocity at LEO (a = 200 km) in km/s

V = sqrt(vinf^2 + v_escape^2);

v_LEO = sqrt(muE/(RE+a));

DV_LEO2EE = abs(V-v_LEO); % km/s
DV_LEO2EE = DV_LEO2EE*10^3; % convert to m/s

DV2 = abs(DV2) * 10^3; % Make sure it's positive and in m/s

m_SATS = computeMaxPayloadMass(maxShuttlePayload, DV_LEO2EE+DV2, Shuttle_dryMass, Shuttle_wetMass-Shuttle_dryMass, Shuttle_Isp, g0);    
Shuttle_FE = []; % todo if still relevant
leftoverShuttle_FM = [];

if ~isnan(m_SATS)
    % Find mass of each satellite
    
    [m_persat, ~, Sat_FEs] = computeSatMForPhasingDeployment(aMOG, eMOG, mu, NSats, Sat_Isp, m_SATS, frac);
    
    Sats_MOGMs = repmat(m_persat, NSats, 1);
    exitflag = 1;
else
    Sats_MOGMs = NaN;
    Sat_FEs = NaN;
    exitflag = 0; % Infeasible launch and insertion
end
end