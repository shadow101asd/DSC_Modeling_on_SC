function [] = DSC7_SC_run(idx, run_idx)

% 1 Circular Loop and 2 Distinct MOGs around E and M State Space Search 
% (very guided by some assumptions): 
% - MOG loops are centered on E, M (fixed center) and each have one
% cartwheel
% Circular loop "bridges" the gap between the two MOGs

assert(idx >= 10, "Index is too small for this architecture to be meaningful!")

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");
[aEai, ~, ~, ~, wEai, ~] = Cartesian2Keplerian(XEa(:,1), muSu);


NSats = SATNUMS(idx)

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% X = [NE eE NM eM aC]

% Constraints!

intcon = [1,3];
lb = [0     0.0  0     1     0.0  0.0  0     1.0];
ub = [NSats 0.6  NSats NSats 0.6  0.6  NSats 1.6];

options = optimoptions('ga', 'Display', 'iter', ...
    'FunctionTolerance', 1e-4, "PopulationSize", 200, ...
    'MaxStallGenerations', 10, 'CreationFcn', 'gacreationuniformint');

nvars = 10;

% Running the GA:

[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFunc7(X,XEa,XMa,etR,muSu,aEai,wEai,aMai,wMai), ...
                                        nvars,[],[],[],[],lb,ub,@(X)nonlcon(X,NSats),intcon,options);

X_opt

% Saving results

filename = "../Data/run"+run_idx+"/" + int2str(NSats);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT");

% Nested Functions

function Out = wrapperFunc7(X,X1,X2,etR,mu,aEai,w0Eai,aMai,w0Mai)
    [NbrE, NshE, eminE, emaxE, NbrM, NshM, eminM, emaxM, NC, aC] = unpackVars(X);

    % Circular Loop
    XSatsC = NSATSpropagateFromKepleriansCIRC(aC,pi,mu,etR,NC);

    % MOGs
    % Earth
    KiE = [aEai; 0.25; 0.0; 0.0; w0Eai; 0.0]; % Second element (eccentricity) isn't actually used here
    shellsE = NbrE*ones(NshE,1);
    esE = linspace(eminE,emaxE,NshE);
    XSatsE = NSATSpropagateFromKepleriansSHELLS(KiE,mu,etR,shellsE,esE,1);
    % Mars
    KiM = [aMai; 0.25; 0.0; 0.0; w0Mai; 0.0]; % Second element (eccentricity) isn't actually used here
    shellsM = NbrM*ones(NshM,1);
    esM = linspace(eminM,emaxM,NshM);
    XSatsM = NSATSpropagateFromKepleriansSHELLS(KiM,mu,etR,shellsM,esM,1);

    % Merge
    XSats = cat(3, XSatsE, XSatsM, XSatsC);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [NbrE, NshE, eminE, emaxE, NbrM, NshM, eminM, emaxM, NC, aC] = unpackVars(X)
    AU = 1.496e8; % 1 AU in km
    NbrE = X(1);
    NshE = X(2);
    eminE = X(3);
    emaxE = X(4);
    NbrM = X(5);
    NshM = X(6);
    eminM = X(7);
    emaxM = X(8);
    NC = X(9);
    aC = X(10) * AU;
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

function [c,ceq] = nonlcon(X,NSats)
    [NbrE, NshE, ~, ~, NbrM, NshM, ~, ~, NC, ~] = unpackVars(X);
    numsats = NbrE*NshE + NbrM*NshM + NC;
    ceq = []; % No equality constraints for MINLPs
    c = [numsats - NSats, NSats - numsats]; % Hacky way of enforcing NSats?
end
end

