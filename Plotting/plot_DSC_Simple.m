%% DSC Simple Plotting

run_idx = "008";
architectures = ["A1", "A2", "A3", "A4", "B", "C1", "C2", "D1", "D2"];
data_start = 10;
data_end = 600;

%% Data Collection

AU = 1.496e8;

input_file = "../Inputs/run" + run_idx + ".mat";
data_file_header = "../Data/run" + run_idx + "/";

load(input_file)
[~,N] = size(SATNUMS);
Na = length(architectures);

Dequivs = NaN([Na, N]);
Dpaths = NaN([Na, N]);

for i = 1:length(architectures)
    a = architectures(i);
    for n = data_start:data_end
        NSATS = SATNUMS(n);
        try 
            filename = data_file_header+a+"/"+int2str(NSATS)+".mat";
            load(filename, "fval", "XSats_opt")
            % [graphslist,numedges,paths,dpaths,~] = networkAnalysis(XEa,XMa,"Earth","Mars",XSats_opt);
    
            Dequivs(i,n) = fval;
            % Dpaths(i,n) = mean(dpaths)/AU;
        catch ME
            warning("No data for architecture " + a + ", N = " + int2str(NSATS))
        end
    end
end

%% Consolidate by type

bestA = min(Dequivs(1:4,:), [], 1);
B = Dequivs(5,:);
bestC = min(Dequivs(6:7,:), [], 1);
bestD = min(Dequivs(8:9,:), [], 1);

%% Analysis

realD = mean(distanceBetweenXs(XEa,XMa)/AU);

% options = fitoptions( 'Method', 'NonlinearLeastSquares',...
%     'Lower',[0 -1e10 0], 'Upper',[1e10 0 1e10]);
% 
% fitA = fit(SATNUMS(10:end)', bestA(10:end)','power2', options);
% fitB = fit(SATNUMS(10:end)', B(10:end)','power2', options);


fitA = fit(SATNUMS(data_start:data_end)', bestA(data_start:data_end)','power2');
fitB = fit(SATNUMS(data_start:data_end)', B(data_start:data_end)','power2');
fitC = fit(SATNUMS(data_start:data_end)', bestC(data_start:data_end)','power2');
fitD = fit(SATNUMS(data_start:data_end)', bestD(data_start:data_end)','power2');

%% Plotting

% Equivalent Distance

% clf(1)

figure(1)
% yline(realD)
hold on
for i = 1:length(architectures)
    a = architectures(i);
    plot(SATNUMS, Dequivs(i,:))
end
hold off
xlabel("Number of relay satellites")
ylabel("Equivalent Distance [AU]")
% legend(["Earth-Mars Average Distance",architectures])
legend(architectures)
title("Performance of Different Simple DSC Architectures for Earth-Mars Downlink")

% ED Consolidated

% clf(2)

figure(2)
plot(SATNUMS, bestA);
hold on
plot(SATNUMS, B);
plot(SATNUMS, bestC);
plot(SATNUMS, bestD);
% plot(SATNUMS(data_start:data_end), fitA(SATNUMS(data_start:data_end)));
% plot(SATNUMS(data_start:data_end), fitB(SATNUMS(data_start:data_end)));
% plot(SATNUMS(data_start:data_end), fitC(SATNUMS(data_start:data_end)));
% plot(SATNUMS(data_start:data_end), fitD(SATNUMS(data_start:data_end)));
hold off
xlabel("Number of relay satellites")
ylabel("Equivalent Distance [AU]")
legend(["Best A", "B", "Best C", "Best D"]);
% legend(["Best A", "B", "Best C", "Best D", "FitA", "FitB", "FitC", "FitD"]);
title("Performance of Different Archetypes of DSC Architectures for Earth-Mars Downlink")

% 
% % Transmission Time
% clf(2)
% 
% c = 3e6/AU*60; % Speed of light (AU/min)
% 
% figure(3)
% yline(realD/c)
% hold on
% for i = 1:length(architectures)
%     a = architectures(i);
%     plot(SATNUMS, movmean(Dpaths(i,:)/c, 10));
% end
% 
% hold off
% xlabel("Number of relay satellites")
% ylabel("Average transmission time [minutes]")
% legend(["Earth-Mars Average Transmission Time",architectures])
% title("Performance of Different Simple DSC Architectures for Earth-Mars Downlink")