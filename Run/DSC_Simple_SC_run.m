function [] = DSC_Simple_SC_run(idx, run_idx, architectures)

% DSC Simple!

assert(idx >= 10, "Index is too small for the simple architectures to be meaningful!")

options = optimoptions('ga', 'Display', 'iter', ...
    'FunctionTolerance', 1e-4, "PopulationSize", 100, ...
    'MaxStallGenerations', 10, 'CreationFcn', 'gacreationuniformint');

for i = 1:length(architectures)
    architecture_str = "DSC_Simple_"+architectures(i)+"_SC_run";
    disp(architecture_str)
    FUNC = str2func(architecture_str);
    FUNC(idx, run_idx, options)
end

