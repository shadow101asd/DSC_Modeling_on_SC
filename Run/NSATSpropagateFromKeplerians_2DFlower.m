function XSats = NSATSpropagateFromKeplerians_2DFlower(Ki,numsats,W,M,etR,mu)
%NSATSPROPAGATEFROMKEPLERIANSCIRC

% if nargin <= 5 
%     % Default assumption is a full formation % To do later...
%     W_frac = 1; 
%     F_frac = 1;
% end

if numsats >= 1
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);
    Kns = KnsfromK_2DFlower(Ki,numsats,mu,W,M);

    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end
    
else
    XSats = []; % Return empty set if no sats are actually requested
end
end
