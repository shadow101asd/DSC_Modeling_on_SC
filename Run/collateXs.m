function Xs = collateXs(X1,X2,XSats)
%COLLATEXS 
    Xs = cat(3, XSats, X1, X2);
end

