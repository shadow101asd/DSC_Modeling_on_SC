function XSats = NSATSpropagateFromKepleriansCIRC(a,f0,mu,etR,numsats)
%NSATSPROPAGATEFROMKEPLERIANSCIRC
if numsats >= 1
    [~, nT] = size(etR);
    XSats(:,:,numsats) = zeros(6,nT);
    Kns = KnsCirc(a,f0,numsats);

    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end
    
else
    XSats = []; % Return empty set if no circs are actually requested
end

    function Kns = KnsCirc(a,f0,numsats)
        df0 = 2*pi/numsats;
        Kns = zeros([6 numsats]);
        for i = 1:numsats
            f0n = f0 + df0*(i-1);
            Kns(:,i) = [a 0 0 0 0 f0n];
        end
    end
end
