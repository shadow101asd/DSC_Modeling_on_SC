function [] = DSC4_SC_run(idx, run_idx)

% Run details
% Circular loops optimization

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

MAXSATS = SATNUMS(idx)
% Num timesteps
[~,nT] = size(XEa);

% Genetic Algorithm

% Constraints!

[aEai, ~, ~, ~, wEai, ~] = Cartesian2Keplerian(XEa(:,1),muSu);
Ki = [aEai; 0.25; 0.0; 0.0; 1.7412; 0.0];

A = [1 -1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
b = [0 0 0 0 0];
intcon = [3,4,5];
lb = [0.0 0.0  1       1       1];
ub = [0.6 0.6  MAXSATS MAXSATS MAXSATS];

options = optimoptions('ga','Display','iter','FunctionTolerance', 1e-3, ...
         'MaxStallGenerations', 7, 'CreationFcn', 'gacreationuniformint',...
         'PopulationSize', 200, 'UseParallel', false);

% Running the GA:
nvars = 5;
    
[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFunc3(XEa,XMa,Ki,X,etR,muSu), ...
                                        nvars,A,b,[],[],lb,ub,@(X)nonlcon(X,MAXSATS),intcon,options);
    
X_opt

% Saving results

filename = "../Data/run"+run_idx+"/" + int2str(MAXSATS);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT")

% Nested Functions

function Out = wrapperFunc3(X1,X2,Ki,X,etR,mu)
    
    [es,shells,cartwheels] = unpackVars(X);
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,shells,es,cartwheels);

    Out = bestLinkBudget(X1,X2,XSats);
end

function [es,shells,cartwheels] = unpackVars(X)
    % Unpack variables
    emin = X(1);
    emax = X(2);
    numshells = X(3);
    numbranches = X(4);
    cartwheels = X(5);
    es = linspace(emin,emax,numshells);
    shells = numbranches*ones(numshells,1);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

function [c,ceq] = nonlcon(X,MAXSATS)
    [~,shells,cartwheels] = unpackVars(X);
    numsats = sum(shells)*cartwheels;
    ceq = []; % No equality constraints for MINLPs
    c = numsats - MAXSATS;
end

end

