function [] = DSC6_SC_run(idx, run_idx)

% 2 Circular Loops and 2 MOG State Space Search (very guided by some
% assumptions): 
% - MOG loop has 1 shell, 2 cartwheels
% MOG loop "bridges" the gap between the two circular loops

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");


NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [Nbranches aMOG eMOG NC1 aC1 aC2]

% Constraints!

A = [0 -1 0 0 1 0; 0 1 0 0 0 -1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
b = [0 0 0 0 0 0];
intcon = [1,4];
lb = [1         0.5  0.0  0.2*NSats 0.5  1.0];
ub = [0.1*NSats 1.5  0.8  0.4*NSats  1.5  1.8];

options = optimoptions('ga', 'Display', 'iter', ...
    'FunctionTolerance', 1e-4, "PopulationSize", 100, ...
    'MaxStallGenerations', 10, 'CreationFcn', 'gacreationuniformint');

nvars = 6;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFunc6(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,options);

X_opt
% Saving results

filename = "../Data/run"+run_idx+"/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT");

% Nested Functions

function Out = wrapperFunc6(X,X1,X2,etR,mu,Ntotal)
    [Nbranches, aMOG, eMOG, NC1, aC1, aC2] = unpackVars(X);
    NC2 = Ntotal - 2*Nbranches - NC1;

    % MOG
    Ki = [aMOG; 0.25; 0.0; 0.0; pi; 0.0]; % Second element (eccentricity) isn't actually used here
    XSatsMOG = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,Nbranches,eMOG,2);

    % Circular Loops
    XSatsC1 = NSATSpropagateFromKepleriansCIRC(aC1,pi,mu,etR,NC1);
    XSatsC2 = NSATSpropagateFromKepleriansCIRC(aC2,pi,mu,etR,NC2);

    % Merge
    XSats = cat(3, XSatsMOG, XSatsC1, XSatsC2);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [Nbranches, aMOG, eMOG, NC1, aC1, aC2] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    Nbranches = X(1);
    aMOG = X(2) * AU;
    eMOG = X(3);
    NC1 = X(4);
    aC1 = X(5) * AU;
    aC2 = X(6) * AU;
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

