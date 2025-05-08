function DVcoeff = MOGPhasingNumericFitNegative(aMOG, eMOG)
% Version of MOGPhasingNumericFitPositive for small negative MOG phase
% changes
    AU = 149600000; % km
    aMOG = aMOG/AU;

    % Poly3 Fit (d = 0)
    a = 10.8756;
    b = -4.5272;
    c = 25.6050;

    DVcoeff = (a*eMOG.^3 + b*eMOG.^2 + c*eMOG) .* aMOG.^(-0.5);
end