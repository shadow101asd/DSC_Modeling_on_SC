function netAs = computeNetPerturbations(XSats, XPs, mus, etR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nsats = size(XSats, 3);
nT = size(etR, 2);
As = zeros([1, nT, nsats]);
nPs = size(XPs, 3);
for i = 1:nsats
    XSat = XSats(:,:,i);
    Ds = distanceBetweenXs(XSat, XPs);
    Avec = zeros([3, nT, 1]);
    for p = 1:nPs
        mu = mus(p);
        dvecs = XSat(1:3,:) - XPs(1:3,:,p);
        dvecs = normalize(dvecs, 1);
        a = mu*dvecs./(Ds(:,:,p).^2);
        Avec = Avec + a;
    end
    Anorm = sqrt(Avec(1,:).^2 + Avec(2,:).^2 + Avec(3,:).^2);
    As(:,:,i) = Anorm;
    % Resulting units: km/s^2.
end

netAs = zeros([nsats, nT]);
for r = 1:nsats
    netAs(r,:) = As(1,:,r);
end
end