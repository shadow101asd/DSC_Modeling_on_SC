function [DV, ToP, tODetails] = computePhasingDVCIRC(R, df, mu)
    % A positive df indicates increasing the true anomaly

    T1 = 2*pi*sqrt(R^3/mu); % Period of circular orbit [s]
    v1 = sqrt(mu/R);

    % Compute characteristics of elliptical transfer orbit
    t = T1*df/(2*pi);
    T2 = T1 - t; % Positive t ensures "speed up" in orbit

    a2 = (sqrt(mu)*T2/(2*pi))^(2/3); % Transfer orbit semimajor axis [km]
    
    if df >= 0
        ra = R;
        rp = 2*a2 - ra;
    elseif df < 0
        rp = R;
        ra = 2*a2 - rp;
    end
    
    h2 = sqrt(2*mu*ra*rp/(ra+rp)); % Transfer orbit angular momentum
    v2 = h2/R;
    e2 = (ra-rp)/(ra+rp);

    % Return computed values
    DV = 2*abs(v1-v2); % x2 for the entry and exit burn of the transfer orbit
    ToP = T2;

    % Optional struct - WARNING: distances in AU, period in days
    AU = 149600000;
    day = 24*3600;
    tODetails = struct;
    tODetails.e = e2;
    tODetails.rp = rp/AU;
    tODetails.ra = ra/AU;
    tODetails.a = a2/AU;
    tODetails.T = T2/day;
end