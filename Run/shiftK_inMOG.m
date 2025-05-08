function K2 = shiftK_inMOG(K1, dPhase, mu)
    K2 = K1;
    [a,e,i,Om,w,f0i] = unpackKeplerian(K2);
    if e > 1e-4 % Orbit non-circular, change argument of periapsis w
        K2(5) = w + dPhase;
    end
    % Regardless, stagger mean anomaly linearly out in time using
    % updateTrueAnomaly
    T = 2*pi*sqrt(a^3/mu); % Orbital Period
    dt = T*(1-dPhase/(2*pi));
    K2(6) = updateTrueAnomaly(a,e,i,Om,w,f0i,mu,dt);
end