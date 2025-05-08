function [] = DSC5_SC_run(idx, run_idx)

% Exhaustive Combinatorics-based Search of the Single-loop of MOGs space

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");
load("../Inputs/MemoTables/3Factors.mat", "ThreeFs");

NSats = SATNUMS(idx);
table_3fs = ThreeFs.("trip"+int2str(NSats));
clear ThreeFs

% Num configs
nC = size(table_3fs, 1);

% Num timesteps
[~,nT] = size(XEa);

% Optimization Algorithms

% Setting up the structure
s = struct();
s.bestfval = Inf;
s.bestConfig = [];
s.bestXopt = [];

% Constraints!

lb = [0.5 0.0  0.0  0];
ub = [1.5 0.6  0.6  2*pi];

options = optimoptions('particleswarm', 'Display', 'iter', ...
    'FunctionTolerance', 1e-4, "SwarmSize", 100, "MaxStallIterations", 10);

nvars = 4;

% Running the PS:

for i = 1:nC
    configNums = table_3fs(i,:);
    config_str = "S" + int2str(configNums(1)) + "B" + int2str(configNums(2)) + "C" + int2str(configNums(3))

    [X_opt, fval, ~, OUTPUT] = particleswarm(@(X) wrapperFunc5(X,XEa,XMa,etR,muSu,configNums),nvars,lb,ub,options);
    
    % Saving this config's performance in the struct
    s.(config_str).Xopt = X_opt;
    s.(config_str).fval = fval;
    s.(config_str).output = OUTPUT;

    % Saving the best for this satnum
    if fval < s.bestfval
        s.bestfval = fval;
        s.bestConfig = configNums;
        s.bestXopt = X_opt;
    end
end

% Saving results

filename = "../Data/run"+run_idx+"/" + int2str(NSats);
X_opt = s.bestXopt;
fval  = s.bestfval;
bestConfig = s.bestConfig;
save(filename, "s", "X_opt", "fval", "bestConfig");

% Nested Functions

function Out = wrapperFunc5(X,X1,X2,etR,mu,configNums)
    AU = 1.496e8; % 1 AU in km
    [es, shells, cartwheels] = unpackConfig(configNums, X);
    Ki = [X(1)*AU; 0.25; 0.0; 0.0; X(4); 0.0]; % Second element (eccentricity) isn't actually used here
    XSats = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,shells,es,cartwheels);
    Out = bestLinkBudget(X1,X2,XSats);
end

function [es, shells, cartwheels] = unpackConfig(configNums, X)
    emin = min(X(2), X(3));
    emax = max(X(2), X(3));
    numshells = configNums(1);
    numbranches = configNums(2);
    cartwheels = configNums(3);
    es = linspace(emin,emax,numshells);
    shells = numbranches*ones(numshells,1);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

end

