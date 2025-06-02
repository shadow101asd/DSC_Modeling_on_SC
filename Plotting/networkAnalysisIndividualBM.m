
% Individual run benchmarking

AU = 149600000; % AU in km

run_idx = '008';
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "etR");

Nsats = 500;
arch  = "B";
load("../Data/run"+run_idx+"/"+arch+"/"+int2str(Nsats)+".mat", "XSats_opt");

%% Testing

[~,~,~,Dpaths_old,Dequivs_old] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats_opt);

[~,~,~,Dpaths_fast,Dequivs_fast] = networkAnalysis_FAST(XEa,XMa,"Earth","Mars",XSats_opt, 1.2);

[~,~,~,Dpaths_BIN,Dequivs_BIN] = networkAnalysis_BIN(XEa,XMa,"Earth","Mars",XSats_opt);
 
[~,~,~,Dpaths_BIN_fast,Dequivs_BIN_fast] = networkAnalysis_BIN_fast(XEa,XMa,"Earth","Mars",XSats_opt);

% [Dequivs, D_bottleneck_avg] = networkAnalysis_Gfast(XEa,XMa,"Earth","Mars",XSats_opt,1.2);

%% Plotting

figure(1)
plot(etR, Dequivs_old/AU)
hold on
plot(etR, Dequivs_fast/AU)
plot(etR, Dequivs_BIN/AU)
plot(etR, Dequivs_BIN_fast/AU)
hold off
legend("Original Approach", "Original Approach with Performance Tweaks", "Binary Search", "Binary Search with Performance Tweaks")