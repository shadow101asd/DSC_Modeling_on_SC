function [] = DSC_Simple_C1_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture C1
% 1 Circular Loop and a Mars MOG State Space Search (very guided by some
% assumptions)

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [NC aC eMOGMa]

% Constraints!

A = [];
b = [];
intcon = 1;
lb = [0      0.8  0  ];
ub = [NSats  1.5  0.4];

nvars = 3;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncC1(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncC1(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/C1/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncC1(X,X1,X2,etR,mu,Ntotal)
    [NC, aC, eMOG2] = unpackVars(X);
    % Compute sat allocations to central ring
    Nb2 = Ntotal - NC;

    % Planet 2 MOG
    XMOG2 = generateMOGNearPlanet(X2, mu, etR, eMOG2, 1, Nb2);

    % Circular Ring
    XSatsC = NSATSpropagateFromKepleriansCIRC(aC,pi,mu,etR,NC);

    % Merge
    XSats = cat(3, XMOG2, XSatsC);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NC, aC, eMOGMa] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    NC = X(1);
    aC = X(2) * AU;
    eMOGMa = X(3);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

