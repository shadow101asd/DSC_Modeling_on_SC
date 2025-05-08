function XSats = generateXsAlongOrbit(XPlanet, mu, etR, numsats)
%GENERATEKNSALONGORBIT
%Similar to NSATpropagateFromKeplerians, but to place satellites along the real orbit of
%planets

if numsats >= 1
    [a_d, enorm_d, i_d, Om_d, w_d, f0_d] = Cartesian2Keplerian(XPlanet, mu);
    Kis = [a_d; enorm_d; i_d; Om_d; w_d; f0_d];
    Ki = mean(Kis, 2);
    [a,e,i,Om,w,~] = unpackKeplerian(Ki);
    f0i = f0_d(1);
    
    Kns(:,numsats) = zeros(6,1);
    T_orbit = 2*pi*sqrt(a^3/mu);
    time_interval = T_orbit/(numsats+1);
    
    for n = 1:numsats
        Kn = Ki;
        
        dt = n*time_interval;
        % Orbit is non-circular, let's stagger the mean anomaly
        Kn(6) = updateTrueAnomaly(a, e, i, Om, w, f0i, mu, dt);
        % Kn(6) = f0i + n*spacing;
        Kns(:,n) = Kn;
    end
    
    % Obtain Xs
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);
    
    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end

else
    XSats = []; % Return empty set if no sats are actually requested
end
