% 15.083 Final Project

AU = 149600000; % AU in km

run_idx = '008';
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "etR");

% numsats = 10:2:500;
numsats = 10:5:300;
N = length(numsats);
architectures = ["A1", "B", "D1"];
A = length(architectures);

% Initialize fields

Dpath_avgs_old = NaN([N,A]);
Dpath_avgs_fast = NaN([N,A]);
Dpath_avgs_BIN = NaN([N,A]);
Dpath_avgs_BIN_fast = NaN([N,A]);

Dequiv_avgs_old = NaN([N,A]);
Dequiv_avgs_fast = NaN([N,A]);
Dequiv_avgs_BIN = NaN([N,A]);
Dequiv_avgs_BIN_fast = NaN([N,A]);

times_old = NaN([N,A]);
times_fast = NaN([N,A]);
times_BIN = NaN([N,A]);
times_BIN_fast = NaN([N,A]);

%% Run with JIT-stabilized benchmarking

for a = 1:A
    arch = architectures(a);

    for i = 1:N
        Nsats = numsats(i);
        load("../Data/run" + run_idx + "/" + arch + "/" + int2str(Nsats) + ".mat", "XSats_opt");

        % JIT WARMUP for each method
        for warm = 1:3
            [~,~,~,~,~] = networkAnalysis(XEa, XMa, "Earth", "Mars", XSats_opt);
            [~,~,~,~,~] = networkAnalysis_FAST(XEa, XMa, "Earth", "Mars", XSats_opt, 1.5);
            [~,~,~,~,~] = networkAnalysis_BIN(XEa, XMa, "Earth", "Mars", XSats_opt);
            [~,~,~,~,~] = networkAnalysis_BIN_fast(XEa, XMa, "Earth", "Mars", XSats_opt);
        end

        % TIMED RUNS
        tic;
        [~,~,~,Dpaths_old,Dequivs_old] = networkAnalysis(XEa, XMa, "Earth", "Mars", XSats_opt);
        times_old(i,a) = toc;

        tic;
        [~,~,~,Dpaths_fast,Dequivs_fast] = networkAnalysis_FAST(XEa, XMa, "Earth", "Mars", XSats_opt, 1.5);
        times_fast(i,a) = toc;

        tic;
        [~,~,~,Dpaths_BIN,Dequivs_BIN] = networkAnalysis_BIN(XEa, XMa, "Earth", "Mars", XSats_opt);
        times_BIN(i,a) = toc;

        tic;
        [~,~,~,Dpaths_BIN_fast,Dequivs_BIN_fast] = networkAnalysis_BIN_fast(XEa, XMa, "Earth", "Mars", XSats_opt);
        times_BIN_fast(i,a) = toc;


        % Store data
        Dpath_avgs_old(i,a) = mean(Dpaths_old);
        Dpath_avgs_fast(i,a) = mean(Dpaths_fast);
        Dpath_avgs_BIN(i,a) = mean(Dpaths_BIN);
        Dpath_avgs_BIN_fast(i,a) = mean(Dpaths_BIN_fast);

        Dequiv_avgs_old(i,a) = mean(Dequivs_old);
        Dequiv_avgs_fast(i,a) = mean(Dequivs_fast);
        Dequiv_avgs_BIN(i,a) = mean(Dequivs_BIN);
        Dequiv_avgs_BIN_fast(i,a) = mean(Dequivs_BIN_fast);

        % Progress display
        disp(string((a-1)*N + i) + "/" + string(A*N) + " iterations complete.");
    end
end



%% Plotting
idx = 1:2;
figure(7)
plot(numsats, mean(times_old_PP(1:N,idx), 2))
hold on
plot(numsats, mean(times_fast_PP(1:N,idx), 2))
plot(numsats, mean(times_BIN_PP(1:N,idx), 2))
plot(numsats, mean(times_BIN_fast_PP(1:N,idx), 2))
hold off
xlabel("Number of satellites")
ylabel("Time [s]")
legend("Original Approach", "Original Approach with Performance Tweaks", "Binary Search", "Binary Search with Performance Tweaks")

figure(8)
plot(numsats, Dequiv_avgs_old(1:N,idx)/AU)
hold on
plot(numsats, Dequiv_avgs_fast(1:N,idx)/AU)
plot(numsats, Dequiv_avgs_BIN(1:N,idx)/AU)
plot(numsats, Dequiv_avgs_BIN_fast(1:N,idx)/AU)
hold off
xlabel("Number of satellites")
ylabel("Bottleneck Distance [AU]")
legend("Original Approach", "Original Approach with Performance Tweaks", "Binary Search", "Binary Search with Performance Tweaks")


figure(9)
plot(numsats, Dpath_avgs_old(1:N,idx)/AU)
hold on
plot(numsats, Dpath_avgs_fast(1:N,idx)/AU)
plot(numsats, Dpath_avgs_BIN(1:N,idx)/AU)
plot(numsats, Dpath_avgs_BIN_fast(1:N,idx)/AU)
hold off
xlabel("Number of satellites")
ylabel("Total Path Distance [AU]")
legend("Original Approach", "Original Approach with Performance Tweaks", "Binary Search", "Binary Search with Performance Tweaks")

%% Fitting

fit_old = fit(numsats', mean(times_old_PP(1:N,idx), 2), "power2")
fit_BIN = fit(numsats', mean(times_BIN_PP(1:N,idx), 2), "power2")
fit_BIN_fast = fit(numsats', mean(times_BIN_fast_PP(1:N,idx), 2), "power2")
fit_FAST = fit(numsats', mean(times_fast_PP(1:N,idx), 2), "power2")
