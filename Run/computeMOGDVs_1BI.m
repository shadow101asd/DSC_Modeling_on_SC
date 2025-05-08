function [SS_DV, MOGSat_DV, HToF, LToF, lambert_details, tODetails] = computeMOGDVs_1BI(aMOG,eMOG,mu)
%computeMOGDVs1 Summary of this function goes here
%   Detailed explanation goes here
%   Units:

% Constants
AU = 1.496e8; % 1 AU in km

if aMOG == AU
    % warning("a = 1AU, computing worst case phasing maneuver over 5 years instead of Hohmann transfer")
    % [SS_DV, HToF, tODetails] = computePhasingDVCIRC(aMOG, dP/5, mu);
    error("Doesn't include a=1AU for now...")
else
    % Hohmann transfer for SS
    [SS_DV, HToF] = computeHohmannTransferCirc2Circ(AU, aMOG, mu);
    tODetails = [];
end

% MOG Insertion 1BI
MOGSat_DV = sqrt(mu/aMOG) * sqrt(2*(1-sqrt(1-eMOG^2)));
LToF = 0; % 1BI
lambert_details = [];
end