function [] = DSC_Simple_A4_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture A4! Slightly different from DSC6
% 1 Earth Ring, 1 Mars Ring, and 2 MOG State Space Search (very guided by 
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

% X = [Nbranches aMOG eMOG]

% Constraints!

A = [];
b = [];
intcon = 1;
lb = [0          1.0  0.0];
ub = [0.5*NSats  1.5  0.6];

nvars = 3;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncA4(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[fval, XSats_opt] = wrapperFuncA4(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/A4/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncA4(X,X1,X2,etR,mu,Ntotal)
    [Nbranches, aMOG, eMOG] = unpackVars(X);
    
    % Planet1 orbital perimeter
    [aP1, eP1, ~, ~, ~, ~] = Cartesian2Keplerian(X1(:,1), mu);
    p1 = elliptical_orbit_perimeter(aP1, eP1);

    % Planet2 orbital perimeter
    [aP2, eP2, ~, ~, ~, ~] = Cartesian2Keplerian(X2(:,1), mu);
    p2 = elliptical_orbit_perimeter(aP2, eP2);
    
    % Compute sat allocations to rings (proportional to perimeter)
    NP1 = (Ntotal - 2*Nbranches)/(1+p2/p1);
    NP1 = round(NP1);
    NP2 = Ntotal - 2*Nbranches - NP1;

    % MOG
    Ki = [aMOG; 0.25; 0.0; 0.0; pi; 0.0]; % Second element (eccentricity) isn't actually used here
    XSatsMOG = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,Nbranches,eMOG,2);

    % Planet 1 Loop
    XSatsP1 = generateXsAlongOrbit(X1, mu, etR, NP1);

    % Planet 2 Loop
    XSatsP2 = generateXsAlongOrbit(X2, mu, etR, NP2);

    % Merge
    XSats = cat(3, XSatsMOG, XSatsP1, XSatsP2);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [Nbranches, aMOG, eMOG] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    Nbranches = X(1);
    aMOG = X(2) * AU;
    eMOG = X(3);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

