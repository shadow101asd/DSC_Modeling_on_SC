function XSats = generate2DP_Flower(aF,eF,w0F,f0F,W,M,etR,mu,Nsats)
    % generate2DP_Flower
    
    Ki = [aF, eF, 0, 0, w0F, f0F];
    XSats = NSATSpropagateFromKeplerians_2DFlower(Ki,Nsats,W,M,etR,mu);
end