function [m_dry_persat, TOD, exitflag] = getSatDryMassesDD_F2D(aFLO,eFLO,NM,Nw,mu,NL,Npl,f,Shuttle_Isp,Sat_Isp,Shuttle_wetMass,Shuttle_dryMass,maxShuttlePayload)
% Units
% Outputs:
% m_dry_persat - [kg]
% TOD          - [Earth years (365.25 days)]
% exitflag     - Boolean

if NM < 0
    error("Flower Constellation Retrograde Deployment is Not Yet Supported.")
end

% Constants
AU = 1.496e8; % 1 AU in km
g0 = 9.81; % m/s^2

% Compute shuttle deployment orbit semi-major axis
RD = aFLO*abs(NM/Nw)^(2/3);

% Shuttle Hohmann transfer to Deployment Orbit
% Cheaper Approach: Integrate first Leg of Hohmann into Earth Escape Burn

% Find desired v_inf
[~, Shuttle_TOF, DV1, DV2] = computeHohmannTransferCirc2Circ(AU, RD, mu); % outputs in km/s
vinf = abs(DV1); % Is this correct independent of the direction?

muE = 3.986004e5; % km3/s2
RE = 6378; % km
a = 200; % LEO orbit altitude (km)
v_escape = sqrt(2*muE/(RE + a)); % Escape velocity at LEO (a = 200 km) in km/s

V = sqrt(vinf^2 + v_escape^2);

v_LEO = sqrt(muE/(RE+a));

DV_LEO2EE = V-v_LEO; % km/s
DV_LEO2EE = DV_LEO2EE*10^3; % convert to m/s

DV2 = abs(DV2) * 10^3; % Make sure it's positive and in m/s

m_SATS = computeMaxPayloadMass(maxShuttlePayload, DV_LEO2EE+DV2, Shuttle_dryMass, Shuttle_wetMass-Shuttle_dryMass, Shuttle_Isp, g0);    

if isnan(m_SATS)
    m_dry_persat = NaN;
    TOD = NaN;
    exitflag = 0; % Infeasible launch and insertion (undefined/negative payload mass)
    return
end
% Otherwise, the launch is feasible, and we proceed.

% Compute satellite insertion DV - direct deployment
resolution = 300;
Sat_DV = compute_DIdeployment_SatDV_2DF(Nw,NM,aFLO,eFLO,resolution);
Sat_DV = Sat_DV * 1e3; % convert km/s to m/s

% Compute final dry mass per satellite
m_sat_wet = m_SATS/Npl; % kg
m_dry_persat = m_sat_wet * exp(-Sat_DV/(Sat_Isp*g0));    

exitflag = 1;

% Add up deployment times

Nsats = NL*Npl; % Total number of satellites being deployed across NL launches
year = 365.25*24*3600; % one year in seconds

T_between_sats = f*Nw/Nsats * 2*pi*sqrt(RD^3/mu); % total wait time in s

TOD_seconds = Shuttle_TOF + T_between_sats * (Nsats-1); % Total deployment time, seconds
TOD = TOD_seconds/year; % convert to Earth years
end