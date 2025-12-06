function XSats = NSATSpropagateFromKepleriansSHELLS(Kis,mu,etR,shells,es,cartwheels)
%NSATSPROPAGATEFROMKEPLERIANS: now with the option to pass in different Kis
if shells(1) >= 1 && cartwheels >= 1
    [~, nT] = size(etR);
    numsats_shell = sum(shells);
    
    if size(Kis,1) == 1
        Kis = repmat(Kis', [1,length(shells)]);
    elseif size(Kis,2) == 1
        Kis = repmat(Kis, [1,length(shells)]);
    end

    XSats(:,:,numsats_shell) = zeros(6,nT);
    Kns = [];
    
    for i = 1:length(shells)
        shellnum = shells(i);
        K = [Kis(1,i); es(i); Kis(3:6,i)];
        Kni = KnsfromK_Inverted(K,shellnum,mu);
        Kns = [Kns, Kni];
    end
    
    if cartwheels > 1
        for k = 2:cartwheels % Add cartwheel symmetricals from the Earth-centric cartwheel
            for n = 1:numsats_shell
                Kns2 = [Kns(1:4,n);Kns(5,n)+2*pi*(k-1)/(cartwheels);Kns(6,n)];
                Kns = [Kns,Kns2];
            end
        end
    end

    for n = 1:(numsats_shell*cartwheels)
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,mu,etR);
    end

else
    XSats = []; % Return empty set if no MOGs are actually requested
end

end

