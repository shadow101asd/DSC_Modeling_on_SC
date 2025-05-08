function [] = DSC_Simple_E_SC_run(idx, run_idx, gaoptions)

% DSC_Simple Architecture E
% 1 Earth-referenced loop of symmetrical MOGs, and 1 MOG near Mars
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

% X = [NbEa eMOGEa eMOGMa]

% Constraints!

A = [];
b = [];
intcon = 1;
lb = [1          0.0  0.0];
ub = [0.4*NSats  0.5  0.4];

nvars = 3;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncE(X, XEa, XMa, etR, muSu, NSats), ...
                                        nvars,A,b,[],[],lb,ub,[],intcon,gaoptions);

X_opt

% Saving results

[~, XSats_opt] = wrapperFuncE(X_opt, XEa, XMa, etR, muSu, NSats);

filename = "../Data/run"+run_idx+"/E/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT", "XSats_opt");

% Nested Functions

function [Out, XSats] = wrapperFuncE(X,X1,X2,etR,mu,Ntotal)
    [Nb1, eMOG1, eMOG2] = unpackVars(X);
    
    % Compute sat allocations
    [a1, ~, ~, ~, ~, ~] = Cartesian2Keplerian(X1(:,1), mu);
    angular_width = 2*computeWoffsetFromMOG(a1, eMOG1, mu);
    desired_Nsyms = floor(2*pi/angular_width);
    max_Nsyms = floor(Ntotal/Nb1);

    Nsyms = min(desired_Nsyms, max_Nsyms); % Make sure we comply with Ntotal
    Nb2 = Ntotal - Nsyms*Nb1;


    % Planet 1 MOG & Symmetricals
    XMOG1 = generateMOGNearPlanet(X1, mu, etR, eMOG1, Nsyms, Nb1);

    % Planet 2 MOG
    XMOG2 = generateMOGNearPlanet(X2, mu, etR, eMOG2, 1, Nb2);

    % Merge
    XSats = cat(3, XMOG1, XMOG2);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NbEa, eMOGEa, eMOGMa] = unpackVars(X)
    NbEa = X(1);
    eMOGEa = X(2);
    eMOGMa = X(3);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

