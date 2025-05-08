function [] = DSC_Simple_B_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture B
% 1 Circular Loops and 2 Planet MOG State Space Search (very guided by some
% assumptions): 
% - 1 MOG near Earth
% - 1 MOG near Mars
% Circular loop bridges the gap between the two MOGs

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [NbEa eMOGEa NbMa eMOGMa aC]

% Constraints!

A = [];
b = [];
intcon = [1,3];
lb = [0          0.0  0          0.0  1.0];
ub = [0.4*NSats  0.5  0.4*NSats  0.4  1.5];

nvars = 5;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncB(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncB(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/B/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncB(X,X1,X2,etR,mu,Ntotal)
    [Nb1, eMOG1, Nb2, eMOG2, aC] = unpackVars(X);
    % Compute sat allocations to central ring
    NC = Ntotal - Nb1 - Nb2;

    % Planet 1 MOG
    XMOG1 = generateMOGNearPlanet(X1, mu, etR, eMOG1, 1, Nb1);

    % Planet 2 MOG
    XMOG2 = generateMOGNearPlanet(X2, mu, etR, eMOG2, 1, Nb2);

    % Circular Ring
    XSatsC = NSATSpropagateFromKepleriansCIRC(aC,pi,mu,etR,NC);

    % Merge
    XSats = cat(3, XMOG1, XMOG2, XSatsC);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NbEa, eMOGEa, NbMa, eMOGMa, aC] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    NbEa = X(1);
    eMOGEa = X(2);
    NbMa = X(3);
    eMOGMa = X(4);
    aC = X(5) * AU;
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

