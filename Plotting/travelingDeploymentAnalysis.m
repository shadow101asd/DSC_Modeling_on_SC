% Traveling Deployment Calcs Optimization
clear all

AU = 149600000; % 1 AU in km
muSu = 1.327124400419393e+11;
day = 3600*24; % 24h in seconds

aMOG = AU;
eMOG = 0.1;
NSATs = 10;

% Precalcs
TMOG = 2*pi*sqrt(aMOG^3/muSu);
maxDT = 0.5*TMOG/day; % in days!!!

%% Running
% 1. Minimizing Shuttle DVs

% Init. variables
ub = (NSATs+1)*TMOG/day;
ts = optimvar('ts', [NSATs, 1], Type='continuous', LowerBound=0, UpperBound=ub);

% Init. optimization problem
optim_prob1 = optimproblem('ObjectiveSense','min');

% Define objective function

optim_prob1.Objective = fcn2optimexpr(@(ts) wrapper_shuttlemin(aMOG, eMOG, muSu, ts, 1),ts);

% Define initial point

ts0 = linspace(0, maxDT*0.5*NSATs, NSATs)';
x0 = struct;
x0.ts = ts0;

% Define constraints

% Bounds on gaps
optim_prob1.Constraints.cons1 = ts(2:end) - ts(1:(end-1)) >= 0;
optim_prob1.Constraints.cons2 = ts(2:end) - ts(1:(end-1)) <= maxDT;

% Solve!
% ga_options = optimoptions('ga','Display', 'iter');
% ga_options.InitialPopulationRange = [0; ub];
% [x_opt1, fval1, ~, output1] = solve(optim_prob, x0, "Solver", "ga", "Options", ga_options);
% [x_opt1, fval1, ~, output1] = solve(optim_prob, "Solver", "ga", "Options", ga_options);
% 
% fmincon_options = optimoptions('fmincon', 'Display', 'iter', "MaxFunctionEvaluations", 1e7);
% [x_opt1, fval1, ~, output1] = solve(optim_prob, x0, "Solver", "fmincon", "Options", fmincon_options);


pattern_options = optimoptions('patternsearch', 'Display', 'iter', "MaxIterations", 1e6);
[x_opt1, fval1, ~, output1] = solve(optim_prob1, x0, "Solver", "patternsearch", "Options", pattern_options);

%% 2. Minimizing Sat. DVs

% Init. variables
ub = (NSATs+1)*TMOG/day;
ts = optimvar('ts', [NSATs, 1], Type='continuous', LowerBound=0, UpperBound=ub);

% Init. optimization problem
optim_prob2 = optimproblem('ObjectiveSense','min');

% Define objective function

optim_prob2.Objective = fcn2optimexpr(@(ts) wrapper_satsmin(aMOG, eMOG, muSu, ts, -1),ts);

% Define initial point

ts0 = linspace(0, maxDT*0.5*NSATs, NSATs)';
x0 = struct;
x0.ts = ts0;


% Define constraints

% Bounds on gaps
optim_prob2.Constraints.cons1 = ts(2:end) - ts(1:(end-1)) >= 0;
optim_prob2.Constraints.cons2 = ts(2:end) - ts(1:(end-1)) <= maxDT;

% Solve!
pattern_options = optimoptions('patternsearch', 'Display', 'iter');
[x_opt2, fval2, ~, output2] = solve(optim_prob2, x0, "Solver", "patternsearch", "Options", pattern_options);


%% 3. Multiobjective!

% Init. variables
ub = (NSATs+1)*TMOG/day;
ts = optimvar('ts', [NSATs, 1], Type='continuous', LowerBound=0, UpperBound=ub);

% Init. optimization problem
optim_prob3 = optimproblem('ObjectiveSense','min');

% Define objective function

% optim_prob3.Objective = fcn2optimexpr(@(ts) wrapper_multiobj(aMOG, eMOG, muSu, ts),ts);
optim_prob3.Objective.first = fcn2optimexpr(@(ts) wrapper_shuttlemin(aMOG, eMOG, muSu, ts, 1),ts);
optim_prob3.Objective.second = fcn2optimexpr(@(ts) wrapper_satsmin(aMOG, eMOG, muSu, ts, 1),ts);

% Define initial point

ts0 = linspace(0, maxDT*0.8*NSATs, NSATs)';
x0 = struct;
x0.ts = ts0;

% Define constraints

% Bounds on gaps
optim_prob3.Constraints.cons1 = ts(2:end) - ts(1:(end-1)) >= 0;
optim_prob3.Constraints.cons2 = ts(2:end) - ts(1:(end-1)) <= maxDT;

% Solve!
% pareto_options = optimoptions('paretosearch', 'Display', 'iter');
% [x_opt3, fval3, ~, output3] = solve(optim_prob3, x0, "Solver", "paretosearch", "Options", pareto_options);

gamultiobj_options = optimoptions('gamultiobj', 'Display', 'iter', 'ParetoFraction', 0.50);
[x_opt3, fval3, ~, output3] = solve(optim_prob3, x0, "Solver", "gamultiobj", "Options", gamultiobj_options);

%% Analysis

%% Plotting

%% Helper funcs

function SUM_SHUTTLEDVS = wrapper_shuttlemin(aMOG, eMOG, muSu, ts, direction)
    day = 3600*24; % 24h in seconds
    [Shuttle_DVs, ~] = travelingDeploymentCalcs(aMOG, eMOG, muSu, ts*day, direction);
    SUM_SHUTTLEDVS = sum(Shuttle_DVs);
end

function SUM_SATSDVS = wrapper_satsmin(aMOG, eMOG, muSu, ts, direction)
    day = 3600*24; % 24h in seconds
    [~, Sats_DVs] = travelingDeploymentCalcs(aMOG, eMOG, muSu, ts*day, direction);
    SUM_SATSDVS = sum(Sats_DVs);
end

function [SUM_SHUTTLEDVS, SUM_SATSDVS] = wrapper_multiobj(aMOG, eMOG, muSu, ts, direction)
    day = 3600*24; % 24h in seconds
    [Shuttle_DVs, Sats_DVs] = travelingDeploymentCalcs(aMOG, eMOG, muSu, ts*day, direction);
    SUM_SHUTTLEDVS = sum(Shuttle_DVs);
    SUM_SATSDVS = sum(Sats_DVs);
end