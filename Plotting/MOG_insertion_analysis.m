% MOG Orbital Insertion Comparisons

%clear all
AU = 149600000; % km
muSu = 1.327124400419393e+11; % km3/s2
aMOG = 1.00*AU;

%% Run
es = 1e-2:1e-2:0.99;
% es = 0.1:0.1:0.5;
% DVs_direct = Inf(size(es));
DVs_lambert = Inf(size(es));
DVs_direct_anal = Inf(size(es));

for i = 1:length(es)
    % DVs_direct(i) = find1BurnTransferFromCirc2MOG(aMOG, es(i), muSu);
    % DVs_lambert(i) = findOptimalTFC2MOGSlot(aMOG, es(i), muSu); OLD PROCEDURE
    DVs_lambert(i) = computeHohmannTransferCirc2Ell(aMOG,aMOG,es(i),muSu);
    DVs_direct_anal(i) = compute1BurnDVTransferCirc2MOG(aMOG, es(i), muSu);
    disp(i+"/"+length(es))
end


%% Plotting

figure(666)
plot(es,DVs_lambert)
hold on
plot(es, DVs_direct)
plot(es, DVs_direct_anal)
hold off
xlabel("MOG Eccentricity");
ylabel("Insertion Total DV from Circular Orbit");
legend("Lambert Problem Search", "1 Burn Direct Insertion", "1 Burn Direct Insertion Analytical Calculation");
title("Comparison of DV Requirements per Satellite Injection")

figure(667)
plot(es, DVs_direct - DVs_direct_anal)
xlabel("MOG Eccentricity");
ylabel("Difference in Insertion Total DV from Circular Orbit");
% Result: just some rounding/numerical approximation differences.
% Conclusion: replace old numeric function with analytical function
