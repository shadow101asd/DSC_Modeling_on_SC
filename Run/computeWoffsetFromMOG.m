function woffset = computeWoffsetFromMOG(aMOG, eMOG, mu)
% Calculates "width" of a given MOG
    numsats = 50; % should be grainy enough for a good estimate
    etR = 1104537669; % etR values don't matter

    % Ki = [aMOG, eMOG, 0, 0, 0, 0];
    % Kns = KnsfromK_Inverted(Ki,numsats,mu);
    % 
    % % Obtain Xs
    % [~, nT] = size(etR);
    % Xs_MOG(:,:,numsats) = zeros(6,nT);
    % 
    % for n = 1:numsats
    %     Kn = Kns(:,n);
    %     Xs_MOG(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    % end
    
    XPlanet = Keplerian2Cartesian(aMOG,0,0,0,pi,pi,mu);
    Xs_MOG = generateMOGNearPlanet(XPlanet, mu, etR, eMOG, 1, numsats, false);

    % Obtain longest distance amongst MOG satellites ("width")
    width = 0;
    XSun = [0;0;0;0;0;0];
    for i = 1:numsats
        for j = i:numsats
            d = distanceBetweenXs(Xs_MOG(:,1,i),Xs_MOG(:,1,j));
            if d > width && ~isnan(d)
                ai = distanceBetweenXs(Xs_MOG(:,1,i), XSun);
                aj = distanceBetweenXs(Xs_MOG(:,1,j), XSun);
                a = max(ai,aj);
                width = d;
            end
        end
    end
    
    if width/(2*a) < 1
        woffset = 2*asin(width/(2*a)) /2; % offset is alpha/2 since we want the offset from the center of the MOG
    else
        error("MOG width calculation error!")
    end
end