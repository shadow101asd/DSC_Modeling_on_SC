function XSats = NSATSpropagateFromKepleriansCIRC(a,f0,mu,etR,numsats,frac)
%NSATSPROPAGATEFROMKEPLERIANSCIRC

if nargin <= 5
    frac = 1; % Default assumption is a full circle, for backwards compatibility
end

if numsats >= 1
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);
    Kns = KnsCirc(a,f0,numsats,frac);

    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end
    
else
    XSats = []; % Return empty set if no circs are actually requested
end

    function Kns = KnsCirc(a,f0,numsats,frac)
        df0 = frac*2*pi/numsats;
        Kns = zeros([6 numsats]);
        for i = 1:numsats
            f0n = f0 + df0*(i-1);
            Kns(:,i) = [a 0 0 0 0 f0n];
        end
    end
end
