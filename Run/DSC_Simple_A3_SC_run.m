function [] = DSC_Simple_A3_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture A3! Slightly different from DSC6
% 1 Circular Loop, 1 Earth Ring, and 2 MOG State Space Search (very guided by 
% some assumptions): 
% - MOG loop has 1 shell, 2 cartwheels
% MOG loop "bridges" the gap between the two other loops

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [Nbranches aMOG eMOG aC2]

% Constraints!

A = [0 1 0 -1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
b = [0 0 0 0];
intcon = 1;
lb = [0          0.5  0.0  1.0];
ub = [0.5*NSats  1.5  0.6  1.8];

nvars = 4;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncA3(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncA3(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/A3/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncA3(X,X1,X2,etR,mu,Ntotal)
    [Nbranches, aMOG, eMOG, aC2] = unpackVars(X);
    pC2 = 2*pi*aC2;

    % Planet2 orbital perimeter
    [aP1, eP1, ~, ~, ~, ~] = Cartesian2Keplerian(X1(:,1), mu);
    p1 = elliptical_orbit_perimeter(aP1, eP1);
    
    % Compute sat allocations to rings (proportional to perimeter)
    NP1 = (Ntotal - 2*Nbranches)/(1+pC2/p1);
    NP1 = round(NP1);
    NC2 = Ntotal - 2*Nbranches - NP1;

    % MOG
    Ki = [aMOG; 0.25; 0.0; 0.0; pi; 0.0]; % Second element (eccentricity) isn't actually used here
    XSatsMOG = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,Nbranches,eMOG,2);

    % Circular Loop
    XSatsC2 = NSATSpropagateFromKepleriansCIRC(aC2,pi,mu,etR,NC2);

    % Planet Loop
    XSatsP1 = generateXsAlongOrbit(X1, mu, etR, NP1);

    % Merge
    XSats = cat(3, XSatsMOG, XSatsC2, XSatsP1);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [Nbranches, aMOG, eMOG, aC2] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    Nbranches = X(1);
    aMOG = X(2) * AU;
    eMOG = X(3);
    aC2 = X(4) * AU;
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

