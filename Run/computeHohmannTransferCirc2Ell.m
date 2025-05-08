function [DV, ToF, DV1, DV2, periapsis, worseDV] = computeHohmannTransferCirc2Ell(R1,a2,ecc2,mu)
%computeHohmannTransferCirc2Circ Summary
%   Returns the DV and ToF of a Hohmann transfer between an initial
%   circular orbit and an arrival elliptical orbit.
%   If R1 = a2, ToF is zero since we consider that there is no
%   need for a Hohmann transfer in that case, only an eccentricity change.
%   If e2 is also zero in that case, then DV is zero (no need for H
%   transfer at all).
%   periapsis is a bool that represents whether the optimal transfer
%   computed and returned is to the elliptical orbit periapsis or apoapsis.
%   
%   Units:
%   Inputs:
%   R1, a2: [km]
%   e2    : []
%   mu    : [km3/s2]
%   Outputs:
%   DV, DV1, DV2: [km/s]
%   ToF:          [s]
%   periapsis:    []
%
%   All returned DVs are magnitudes (taken as positive)

if R1 == a2 && ecc2 == 0 % There is no need for a Hohmann Transfer!
    DV = 0.0;
    ToF = 0.0;
    DV1 = 0.0;
    DV2 = 0.0;
    periapsis = false;
elseif ecc2 == 0
    [DV, ToF, DV1, DV2] = computeHohmannTransferCirc2Circ(R1,a2,mu);
    periapsis = false;
else
    % Two cases to check: transfer to periapsis, and transfer to apoapsis
    
    % Periapsis
    R2p = a2*(1+ecc2); % Periapsis distance from Sun [km]
    a_tp = (R1+R2p)/2;
    DV1_p = abs(sqrt(mu/R1)*(sqrt(2*R2p/(R1+R2p)) - 1));

    % From vis viva:
    v_tatrp = sqrt(mu*(2/R2p - 1/a_tp));
    v2p = sqrt(mu*(2/R2p - 1/a2));

    DV2_p = abs(v_tatrp - v2p);
    DV_p = DV1_p + DV2_p;
    
    % Apoapsis
    R2a = a2*(1-ecc2); % Apoapsis distance from Sun [km]
    a_ta = (R1+R2a)/2;
    DV1_a = abs(sqrt(mu/R1)*(sqrt(2*R2a/(R1+R2a)) - 1));
    
    % From vis viva:
    v_tatra = sqrt(mu*(2/R2a - 1/a_ta));
    v2a = sqrt(mu*(2/R2a - 1/a2));

    DV2_a = abs(v_tatra - v2a);
    DV_a = DV1_a + DV2_a;

    % Take cheapest transfer and compute outputs:

    if DV_a <= DV_p
        DV = DV_a;
        DV1 = DV1_a;
        DV2 = DV2_a;
        periapsis = false;
        a_t = a_ta;
        worseDV = DV_p; % might not really be used
    else
        DV = DV_p;
        DV1 = DV1_p;
        DV2 = DV2_p;
        periapsis = true;
        a_t = a_tp;
        worseDV = DV_a; % might not really be used
    end

    % Time of flight. 
    % Transfer orbit orbital period
    T_transfer = 2*pi*sqrt(a_t^3 / mu);
    
    % HT covers half of the elliptical transfer orbit, of
    % period (from perigee to apogee or vice versa)
    ToF = T_transfer/2;
end
end