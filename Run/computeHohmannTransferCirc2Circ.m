function [DV, ToF, DV1, DV2] = computeHohmannTransferCirc2Circ(R1,R2,mu)
%computeHohmannTransferCirc2Circ Summary
%   Returns the DV and ToF of a Hohmann transfer between 2 circular orbits.
%   If R1 = R2, Dvs and ToF are zero since we consider that there is no
%   need for a Hohmann transfer in that case.
%   
%   Units:
%   Inputs:
%   R1, R2: [km]
%   mu    : [km3/s2]
%   Outputs:
%   DV, DV1, DV2: [km/s]
%   ToF:          [s]

if R1 == R2 % There is no need for a Hohmann Transfer!
    DV = 0.0;
    ToF = 0.0;
    DV1 = 0.0;
    DV2 = 0.0;
else
    a_t = (R1+R2)/2;
    DV1 = abs(sqrt(mu/R1)*(sqrt(2*R2/(R1+R2)) - 1));
    DV2 = abs(sqrt(mu/R2)*(1 - sqrt(2*R1/(R1+R2))));

    DV = DV1 + DV2;

    % Time of flight. 
    % Transfer orbit orbital period
    T_transfer = 2*pi*sqrt(a_t^3 / mu);
    
    % HT covers half of the elliptical transfer orbit, of
    % period (from perigee to apogee or vice versa)
    ToF = T_transfer/2;
end
end