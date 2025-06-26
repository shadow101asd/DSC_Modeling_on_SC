function [DVs, ToPs] = getDVsForPhasingDeployment(aMOG, eMOG, mu, NSats, frac)

    if nargin <= 4
        frac = 1; % The default assumption is that we're deploying a full MOG
    end

    spacing = 2*pi/NSats;
    if rem(NSats, 2) == 0 % If NSats is even
        dPs = linspace(spacing/2, 2*pi-spacing/2, NSats);
    else % If NSats is odd
        dPs = linspace(0, 2*pi - spacing, NSats);
    end

    % Scale by frac for partial deployment
    dPs = dPs * frac;
    
    DVs = Inf(size(dPs));
    ToPs = DVs;
    
    for i = 1:NSats
        dP = dPs(i);
        [DV, ToP] = computeOptimalPhasingDVMOG_extended(aMOG, eMOG, dP, mu); % Changed to extended since this is more reliable and faster
        DVs(i) = DV;
        ToPs(i) = ToP;
    end
end