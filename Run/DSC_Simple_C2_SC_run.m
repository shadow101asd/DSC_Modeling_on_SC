function [] = DSC_Simple_C2_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture C2
% 1 Fixed loop at Earth and a Mars MOG State Space Search

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [NEa eMOGMa]

% Constraints!

A = [];
b = [];
intcon = 1;
lb = [0      0];
ub = [NSats  0.4];

nvars = 2;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncC2(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncC2(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/C2/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncC2(X,X1,X2,etR,mu,Ntotal)
    [NP1, eMOG2] = unpackVars(X);
    % Compute sat allocations to central ring
    Nb2 = Ntotal - NP1;

    % Planet 1 Ring
    XSatsP1 = generateXsAlongOrbit(X1, mu, etR, NP1);

    % Planet 2 MOG
    XMOG2 = generateMOGNearPlanet(X2, mu, etR, eMOG2, 1, Nb2);

    % Merge
    XSats = cat(3, XMOG2, XSatsP1);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NEa, eMOGMa] = unpackVars(X)
    NEa = X(1);
    eMOGMa = X(2);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

