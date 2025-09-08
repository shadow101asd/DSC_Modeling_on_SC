function Kns = KnsfromK_2DFlower(Ki,numsats,mu,W,M)
%KNSFROMK but for 2D Flower Constellations

    [a,e,i,Om,w,f0i] = unpackKeplerian(Ki);
    Kns(:,numsats) = zeros(6,1);
    % Assume orbital symmetry around the Sun for now
    for n = 1:numsats
        Kn = Ki;
        spacing_W = 2*pi/numsats * W;
        if e > 1e-4 % Orbit non-circular, change argument of periapsis w and eccentricity if applicable
            Kn(5) = w + spacing_W*(n-1);
        end

        % Regardless, stagger mean anomaly linearly out in time using
        % updateTrueAnomaly
        T = 2*pi*sqrt(a^3/mu); % Orbital Period
        dt = T + M*T*(n-1)/numsats;
        Kn(6) = updateTrueAnomaly(a,e,i,Om,w,f0i,mu,dt);

        Kns(:,n) = Kn;
    end
end