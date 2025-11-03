function K2 = shiftK_in_2DF(K1, W, M, dPhase, mu)
    K2 = K1;
    [a,e,i,Om,w,f0i] = unpackKeplerian(K2);
    if e > 1e-4 % Orbit non-circular, change argument of periapsis w
        K2(5) = w + dPhase * W;
    end
    % Regardless, update mean anomaly using updateTrueAnomaly
    T = 2*pi*sqrt(a^3/mu); % Orbital Period
    dt = T*(1 + M*dPhase/(2*pi));
    K2(6) = updateTrueAnomaly(a,e,i,Om,w,f0i,mu,dt);
end