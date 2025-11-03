function XSats = NSATSpropagateFromKeplerians_2DFlower_Es(Ki,numsats,W,M,E,etR,mu,emin,emax,Es_distribution_name)
%NSATSpropagateFromKeplerians_2DFlower_Es

% if nargin <= 5 
%     % Default assumption is a full formation % To do later...
%     W_frac = 1; 
%     F_frac = 1;
%     E_frac = 1;
% end

if numsats >= 1
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);

    if nargin <= 8
        Kns = KnsfromK_2DFlower_Es(Ki,numsats,mu,W,M,E);
    elseif nargin == 9
        Kns = KnsfromK_2DFlower_Es(Ki,numsats,mu,W,M,E,emin,emax);
    else
        Kns = KnsfromK_2DFlower_Es(Ki,numsats,mu,W,M,E,emin,emax,Es_distribution_name);
    end

    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end
    
else
    XSats = []; % Return empty set if no sats are actually requested
end
end
