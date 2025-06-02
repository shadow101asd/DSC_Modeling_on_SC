function [Sats_MOGMs, Shuttle_FE, Sat_FEs, leftoverShuttle_FM, exitflag] = getSatDryMassesSI(aMOG,eMOG,mu,NSats,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload)
%getSatDryMassesSI Summary of this function goes here
%   Detailed explanation goes here

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

m0 = Shuttle_wetMass; % Is this correct? I think this is supposed to include the payload mass... The current approach is only slightly wrong as far as I can tell... Would have to solve an equation for payload mass here.
mEE = m0  * exp(-DV_LEO2EE/(Shuttle_Isp*g0));

% Find payload mass after second leg of Hohmann transfer

DV2 = abs(DV2) * 10^3; % Make sure it's positive and in m/s
mHA = mEE  * exp(-DV2/(Shuttle_Isp*g0));
Shuttle_FE = m0 - mHA;

% End of new cheaper approach

% Starship Block III can carry up to 200t paylaod to LEO before refuelling
m_SATS = min(maxShuttlePayload, mHA-Shuttle_dryMass); % kg
leftoverShuttle_FM = mHA-Shuttle_dryMass - m_SATS;

if m_SATS > 0
    % Find mass of each satellite
    
    [m_persat, ~, Sat_FEs] = computeSatMForPhasingDeployment(aMOG, eMOG, mu, NSats, Sat_Isp, m_SATS);
    
    Sats_MOGMs = repmat(m_persat, NSats, 1);
    exitflag = 1;
else
    Sats_MOGMs = NaN;
    Sat_FEs = NaN;
    exitflag = 0; % Infeasible launch and insertion (negative payload mass)
end
end