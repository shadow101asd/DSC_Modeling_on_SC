function [D_bottleneck, exitflag] = bottleneckshortestpath_GUROBI(A)
%BOTTLENECKSHORTESTPATH using Gurobi. It is assumed that the source and
%target nodes are the last two indices of the adjacency matrix.

Nsats = size(A, 1) - 2;
s = Nsats + 1; % Source idx
ta = Nsats + 2; % Target idx
N = ta;
A_bin = double(A > 0);

% Construct Optimization Problem

% Init. variables

% 2D Integer Matrix reprensenting whether edge between nodes (i,j) is
% selected in the BSP
X_ij = optimvar('X_ij', [N, N], Type='integer', LowerBound=0, UpperBound=A_bin); 

dmax = optimvar('dmax', Type='continuous', LowerBound=0, UpperBound=inf);

% Init. optimization problem
optim_problem = optimproblem('ObjectiveSense','min');

% Define objective function
optim_problem.Objective = dmax;

% Define constraints

% dmax coupling constraints
optim_problem.Constraints.cons_DMAX = A .* X_ij <= dmax;

% Flow conservation

optim_problem.Constraints.cons_target = sum(X_ij(:,ta), 1) - sum(X_ij(ta,:), 2) ==  1;
optim_problem.Constraints.cons_source = sum(X_ij(:,s), 1) - sum(X_ij(s,:), 2) == -1;

optim_problem.Constraints.cons_sats = sum(X_ij(:,1:Nsats), 1) - sum(X_ij(1:Nsats,:), 2)' == 0;

% Tightening

optim_problem.Constraints.cons_tightening = sum(X_ij, 1) <= 1;   % no node is visited more than once
optim_problem.Constraints.cons_tightening = sum(X_ij, 2) <= 1;   % no node is visited more than once


[~, D_bottleneck, exitflag, output] = solve(optim_problem);
    
end

