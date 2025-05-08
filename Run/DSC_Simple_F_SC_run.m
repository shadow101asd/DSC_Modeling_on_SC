function [] = DSC_Simple_F_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture F
% 1 Mars-referenced loop of symmetrical MOGs, and 1 MOG near Earth
% All MOGs have depth 1, i.e. only one layer

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [eMOGEa NbMa eMOGMa]

% Constraints!

A = [];
b = [];
intcon = 2;
lb = [0.0  1         0.0];
ub = [0.5  0.4*NSats 0.4];

nvars = 3;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncF(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncF(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/F/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncF(X,X1,X2,etR,mu,Ntotal)
    [eMOG1, Nb2, eMOG2] = unpackVars(X);
    
    % Compute sat allocations
    [a2, ~, ~, ~, ~, ~] = Cartesian2Keplerian(X2(:,1), mu);
    angular_width = 2*computeWoffsetFromMOG(a2, eMOG2, mu);
    desired_Nsyms = floor(2*pi/angular_width);
    max_Nsyms = floor(Ntotal/Nb2);

    Nsyms = min(desired_Nsyms, max_Nsyms); % Make sure we comply with Ntotal
    Nb1 = Ntotal - Nsyms*Nb2;


    % Planet 1 MOG
    XMOG1 = generateMOGNearPlanet(X1, mu, etR, eMOG1, 1, Nb1);

    % Planet 2 MOG & Symmetricals
    XMOG2 = generateMOGNearPlanet(X2, mu, etR, eMOG2, Nsyms, Nb2);

    % Merge
    XSats = cat(3, XMOG1, XMOG2);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [eMOGEa, NbMa, eMOGMa] = unpackVars(X)
    eMOGEa = X(1);
    NbMa = X(2);
    eMOGMa = X(3);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

