function [theta, rho] = getPolarRepresentationOfMOG(aMOG, eMOG, mu, Nsteps)
    AU = 1.496e8; % 1 AU in km

    % Optional resolution argument
    if nargin == 3
        Nsteps = 1e4;
    end
    
    etR = 1;
    XP = Keplerian2Cartesian(aMOG, 0, 0, 0, 0, pi/2, mu);
    XSlots = generateMOGNearPlanet(XP, mu, etR, eMOG, 1, Nsteps,false);

    % Shift MOG to be centered around origin
    XSlots = XSlots - XP;
    XSlots = XSlots/AU; % Convert to AU and AU/s
    [theta, rho] = cart2pol(XSlots(1,:,:), XSlots(2,:,:));
    theta = reshape(theta, [Nsteps,1,1]);
    rho = reshape(rho, [Nsteps,1,1]);
end