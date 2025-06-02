function [Dequivs, D_bottleneck_avg, output] = networkAnalysis_G(X1,X2,~,~,XSats)
    AU = 149600000; % AU in km

    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats)/AU;
    D12 = distanceBetweenXs(X1,X2)/AU;
    [~,T,Nsats] = size(XSats);
    A = createAdjacencyMatrix2(Xs,@distanceBetweenXs,D12);

    A_bin = double(A > 0);
    si = Nsats + 1; % Source idx
    ta = Nsats + 2; % Target idx
    N = ta;


    % Construct Optimization Problem

    % Init. variables

    % 2D Integer Matrix reprensenting whether edge between nodes (i,j) is
    % selected in the BSP
    X_ijt = optimvar('X_ijt', [N, N, T], Type='integer', LowerBound=0, UpperBound=A_bin); 

    dmax_t = optimvar('dmax_t', [1,1,T], Type='continuous', LowerBound=0, UpperBound=inf);

    % Init. optimization problem
    optim_problem = optimproblem('ObjectiveSense','min');

    % Define objective function
    optim_problem.Objective = sum(dmax_t)/T; % Average Dequiv
    
    % Define constraints

    % dmax coupling constraints
    optim_problem.Constraints.cons_DMAX = A .* X_ijt <= repmat(dmax_t, N, N, 1);

    % Flow conservation

    optim_problem.Constraints.cons_target = sum(X_ijt(:,ta,:), 1) - sum(X_ijt(ta,:,:), 2) ==  1;
    optim_problem.Constraints.cons_source = sum(X_ijt(:,si,:), 1) - sum(X_ijt(si,:,:), 2) == -1;

    for k = 1:Nsats
        optim_problem.Constraints.(['cons_', num2str(k)]) = sum(X_ijt(:,k,:), 1) - sum(X_ijt(k,:,:), 2) == 0;
    end
    
    % Solve using Gurobi
    % % Set up paths: This needs to be fixed to run on SC or elsehwere than this Mac
    % gurobi_file_path = horzcat('/Library/gurobi1201/macos_universal2', '/matlab/gurobi_setup.m');
    % gurobi_add_path = horzcat('/Library/gurobi1201/macos_universal2', '/examples/matlab/');
    % 
    % % Setup gurobi
    % run(gurobi_file_path)
    % addpath(gurobi_add_path); % temporarily overwrites MATLAB's 'intlinprog.m'
    
    opts = optimoptions('intlinprog', 'Display', 'iter');
    [sol, D_bottleneck_avg, ~, output] = solve(optim_problem, Options=opts);
    D_bottleneck_avg = D_bottleneck_avg * AU;

    Dequivs = sol.dmax_t * AU;
end