function [Sat_DV, theta_intersection, dmin, exitflag] = compute_DIdeployment_SatDV_2DF(Nw,NM,aF,eF,resolution,mu)
    % Returned units: km/s

    if nargin < 5
        resolution = 200;
    end

    if nargin < 6
        mu = 1.327124400419393e+11;
    end

    RD = aF * abs(NM/Nw)^(2/3);
    
    % Return NaN if there is no intersection between the orbits
    if (RD > aF*(1+eF)) || (RD < aF*(1-eF))
        Sat_DV = NaN;
        theta_intersection = NaN;
        dmin = NaN;
        exitflag = 0;
        return
    end

    Kcirc = [RD, 0, 0, 0, 0, 0]';
    Kfc = [aF, eF, 0, 0, 0, 0]';
    
    thetas = linspace(pi, 2*pi, resolution);

    Kcircs = repmat(Kcirc, [1, resolution]);
    Kcircs(6,:) = thetas;

    Kfcs = repmat(Kfc, [1, resolution]);
    Kfcs(6,:) = thetas;
    
    % Convert to Cartesian coordinates
    Xcircs = Keplerian2Cartesian(Kcircs(1,:),Kcircs(2,:),Kcircs(3,:),Kcircs(4,:),Kcircs(5,:),Kcircs(6,:),mu);
    Xfcs = Keplerian2Cartesian(Kfcs(1,:),Kfcs(2,:),Kfcs(3,:),Kfcs(4,:),Kfcs(5,:),Kfcs(6,:),mu);

    % The intersection must occur at the same true anomaly for both!
    [dmin, idx] = min(sqrt(sum((Xcircs(1:3,:)-Xfcs(1:3,:)).^2, 1)));
    theta_intersection = thetas(idx) - pi;
    
    V1 = Xcircs(4:6, idx);
    V2 = Xfcs(4:6, idx);
    Sat_DV = norm(V1-V2, 2);
    exitflag = 1;
end