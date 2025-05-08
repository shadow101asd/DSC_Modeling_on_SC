function [MOGSat_DV, V2, MOGSat_DV_worse, V2_worse] = find1BurnTransferFromCirc2MOG(aMOG,eMOG,mu)
    % OBSOLETE: Use compute1BurnDVTransferCirc2MOG instead where possible

    % Meta parameters
    TMOG = 2*pi*sqrt(aMOG^3/mu); % MOG orbital period
    nT = 1e5;
    etR = linspace(0,2*TMOG/2,nT); % Only half a period needed thanks ot symmetry in the problem

    % Generate Orbital Points
    Ki = [aMOG; eMOG; 0; 0; -pi/2; 0];
    XE = NSATSpropagateFromKepleriansSHELLS(Ki, mu, etR, 1, eMOG, 1);

    % Find intersections of MOG orbit and Circular Orbit
    XSun = zeros(size(XE(:,1)));
    [minError, minIdx] = min(abs(distanceBetweenXs(XSun, XE)-aMOG));
    
    % Get equivalent circ velocities at the intersection
    VC = getCircularOrbitalVelocityAtPos(XE(1:3, minIdx), mu);

    VE = XE(4:6, minIdx);

    % DVs
    MOGSat_DV = abs(norm(VE-VC));
    V2 = VE;
    
    % The symmetry of the problem means we don't need to check both
    % intersections...

    function V = getCircularOrbitalVelocityAtPos(X, mu)
        % Always assuming CCW motion, and equatorial
        X = X(1:3);
        assert(X(3) == 0.0) % equatorial assumed

        R = distanceBetweenXs(X, [0;0;0]);
        s_C = sqrt(mu/R); % km/s
        alpha = atan2(X(2), X(1));
        V = s_C*[-sin(alpha); cos(alpha); 0];
    end
end