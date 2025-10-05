function Xs = propagateFromCartesians(Xis,mu,etR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    nT = length(etR);
    N = size(Xis, 3);
    Xs = zeros([6, nT, N]);
    
    for n = 1:N
        [a, ecc, i, Om, w, f0] = Cartesian2Keplerian(Xis(:,1,n),mu);
        K = [a, ecc, i, Om, w, f0];
        XSat = propagateFromKeplerians(K,mu,etR);
        Xs(:,:,n) = XSat;
    end
end