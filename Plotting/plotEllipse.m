function [] = plotEllipse(K, mu, TMOG_fraction, scaling_factor, colorstring)
    TMOG = 2*pi*sqrt(K(1)^3 /mu);
    etR = linspace(0, TMOG_fraction*TMOG, 1000);
    X = propagateFromKeplerians(K, mu, etR);
    plot(X(1,:)/scaling_factor, X(2,:)/scaling_factor, Color=colorstring)
end