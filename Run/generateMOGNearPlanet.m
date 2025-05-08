function XSats = generateMOGNearPlanet(XPlanet, mu, etR, eMOG, Ncartwheels, numsatsperMOG, offset_bool)
%generateMOGNearPlanet
%Similar to KnsfromK, but to place satellites in a NSats MOG of
%eccentricity eMOG near the planet XPlanet
% No inclination included at the moment

 if ~exist('offset_bool','var')
     offset_bool = true; % we include offset by default
 end


% eMOG of 0 is converted to a small value for numerical stability/continuity
eMOG = max(eMOG, 1e-3);

if numsatsperMOG >= 1 && Ncartwheels >= 1
    [a_d, enorm_d, i_d, Om_d, w_d, f0_d] = Cartesian2Keplerian(XPlanet, mu);
    Kis = [a_d; enorm_d; i_d; Om_d; w_d; f0_d];
    Ki = mean(Kis, 2);
    [a,~,~,Om,w,~] = unpackKeplerian(Ki);
    f0i = f0_d(1);
    
    if offset_bool
        % Compute woffset from provided eMOG and aMOG
        woffset = computeWoffsetFromMOG(a, eMOG, mu);
    else
        woffset = 0.0;
    end
    
    if ~isnan(w) && ~isnan(Om)
        Ki = [a; eMOG; 0; 0; w+Om+f0i+1.0*woffset; 0];
    else
        % Then we're dealing with an idealized orbit
        Ki = [a; eMOG; 0; 0; f0i+1.0*woffset; 0];
    end
    
    % Obtain Xs
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki, mu, etR, numsatsperMOG, eMOG, Ncartwheels);

else
    XSats = []; % Return empty set if no sats are actually requested
end

