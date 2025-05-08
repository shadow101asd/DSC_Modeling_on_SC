function MOGSat_DV = compute1BurnDVTransferCirc2MOG(aMOG,eMOG,mu)
    % See AAS/AIAA MOG Paper for details
    MOGSat_DV = sqrt(mu./aMOG) .* sqrt(2*(1-sqrt(1-eMOG.^2)));
end