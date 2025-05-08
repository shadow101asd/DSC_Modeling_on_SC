function [] = DSC_SC_runN(idx_start, N, period, run_idx, DSC_num)
%DSC4_SC_RUNN Wrapper for Batching N times DSC4_SC_run()
% Runs largest sims first to check memory usage
    FUNC = str2func("DSC"+int2str(DSC_num)+"_SC_run");
    for i = 1:N
        new_idx = idx_start + (N-i)*period;
        FUNC(new_idx, run_idx)
    end
end

