function [DVs, ToPs] = getDVsForPhasingDeployment(aMOG, eMOG, mu, NSats)
    spacing = 2*pi/NSats;
    if rem(NSats, 2) == 0 % If NSats is even
        dPs = linspace(spacing/2, 2*pi-spacing/2, NSats);
    else % If NSats is odd
        dPs = linspace(0, 2*pi - spacing, NSats);
    end
    
    DVs = Inf(size(dPs));
    ToPs = DVs;
    
    for i = 1:NSats
        dP = dPs(i);
        [DV, ToP] = computeOptimalPhasingDVMOG(aMOG, eMOG, dP, mu);
        DVs(i) = DV;
        ToPs(i) = ToP;
    end
end