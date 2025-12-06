function [Xs] = propagateFromKeplerians_plural(Kis,mu,etR)
    
    N = size(Kis,2);
    nT = length(etR);
    Xs = zeros([6,nT,N]);

    for n = 1:N
        Xs(:,:,n) = propagateFromKeplerians(Kis(:,n),mu,etR);
    end

end