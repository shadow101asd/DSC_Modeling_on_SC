function [] = DSC_Simple_D2_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture D2
% 1 Fixed loop at Mars and an Earth MOG State Space Search
% Very similar to C2

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [NMa eMOGEa]

% Constraints!

A = [];
b = [];
intcon = 1;
lb = [0      0];
ub = [NSats  0.4];

nvars = 2;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncD2(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncD2(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/D2/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncD2(X,X1,X2,etR,mu,Ntotal)
    [NP2, eMOG1] = unpackVars(X);
    % Compute sat allocations to central ring
    Nb1 = Ntotal - NP2;

    % Planet 1 MOG
    XMOG1 = generateMOGNearPlanet(X1, mu, etR, eMOG1, 1, Nb1);

    % Planet 2 Ring
    XSatsP2 = generateXsAlongOrbit(X2, mu, etR, NP2);

    % Merge
    XSats = cat(3, XMOG1, XSatsP2);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NMa, eMOGEa] = unpackVars(X)
    NMa = X(1);
    eMOGEa = X(2);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

