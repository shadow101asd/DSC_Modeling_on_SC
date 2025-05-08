function DVcoeff = MOGPhasingNumericFitPositive(aMOG, eMOG)
    AU = 149600000; % km
    aMOG = aMOG/AU;

    % Poly3 Fit (d = 0)
    a = 15.1323;
    b = -4.4752;
    c = 25.8481;

    DVcoeff = (a*eMOG.^3 + b*eMOG.^2 + c*eMOG) .* aMOG.^(-0.5);
end