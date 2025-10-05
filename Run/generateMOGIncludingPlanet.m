function XSats = generateMOGIncludingPlanet(XPlanet, mu, etR, numsatsperMOG)
%generateMOGNearPlanet
%Similar to KnsfromK, but to place satellites in a NSats + 1 MOG that
%includes the planet described by XPlanet

% Should be robust to both real ephemeride data and idealized planar cases.

if numsatsperMOG >= 1
    [aP, eP, iP, OmP, wP, f0P] = Cartesian2Keplerian(XPlanet(:,1), mu);

    if ~isnan(wP) && ~isnan(OmP)
        KPi = [aP; eP; iP; OmP; wP; f0P];
    else
        % Then we're dealing with an idealized orbit
        KPi = [aP; eP; 0; 0; f0P; 0];
    end
    
    XSats_P = NSATSpropagateFromKepleriansSHELLS(KPi, mu, etR, numsatsperMOG+1, eP, 1);
    XSats = XSats_P(:,:,2:end);

else
    XSats = []; % Return empty set if no sats are actually requested
end

