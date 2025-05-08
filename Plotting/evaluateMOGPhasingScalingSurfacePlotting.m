%% Running


AU = 149600000; % km
% as = AU*(0.7:0.05:1.5);
% es = 1e-2:1e-2:5e-1;
as = AU*1.0;
es = linspace(1-1e-5, 1-1e-15, 2);

coeffs = NaN([length(as), length(es)]);
for i = 1:length(as)
    for j = 1:length(es)
        c = linearfit(as(i), es(j), -1);
        coeffs(i,j) = abs(c);
        disp("Progress: " + num2str((i-1)*length(as) + j) + "/" + num2str(length(as)*length(es)) + " analyses")
    end
end

%% Analysis
% options = fitoptions( 'Method', 'NonlinearLeastSquares',...
%     'Lower',[0], 'Upper',[100], 'StartPoint', 5);
% fit1 = fit(dPs(2:end)', DVs(2:end)', 'a*x')

% fit1 = fit(dPs(1:end)', DVs(1:end)', 'poly2')

%% Plotting

% figure(6666)
% plot(dPs, DVs)
% hold on
% plot(dPs, fit1(dPs))
% hold off
% xlabel("MOG Phase Shift [rad]");
% ylabel("Optimal DV Found [km/s]");
% legend("Data", "Fit")
% 
% day = 24*3600;
% figure(6667)
% plot(dPs, ToPs/day)
% xlabel("MOG Phase Shift [rad]");
% ylabel("Time of Flight for Phase Shift (Earth days)");

% figure(6668)
% plot(es, coeffs)
% xlabel("MOG eccentricity");
% ylabel("Slope of resulting distribution [km/s/rad]");

% [X,Y] = meshgrid(es, as/AU);
% figure(6669)
% surf(X,Y,coeffs)
% xlabel("MOG eccentricity");
% ylabel("MOG semi-major axis [AU]")
% zlabel("Cost to make small phase change [km/s/rad]");

figure(6670)
clf(6670)
LEGEND = [];
hold on
plot(es, coeffs(1:4:length(as),:), "-o")
for i = 1:4:length(as)
    LEGEND = [LEGEND, "aMOG = "+ num2str(as(i)/AU) + " AU"];
end
xlabel("MOG eccentricity");
ylabel("Slope of resulting distribution [km/s/rad]");
legend(LEGEND)
hold off

% figure(6671)
% hold on
% for i = 1:length(as)
%     plot(as/AU, coeffs)
% end
% xlabel("aMOG [AU]");
% ylabel("Slope of resulting distribution [km/s/rad]");
% hold off

%% Fitting

fits = struct;
gofs = struct;
for i = 1:length(as)
    %[fiti, gofi] = fit(es', coeffs(i,:)', "a*x/(1-x^b)");
    [fiti, gofi] = fit(es', coeffs(i,:)', "a*x^3 + b*x^2 + c*x");
    fits.("fit"+i) = fiti;
    gofs.("gof"+i) = gofi;
end

figure(27497234)
scatter(es', coeffs(i,:)')
hold on
plot(es', fiti(es'))
% plot(0.01:0.01:0.99, fit1(0.01:0.01:0.99))
hold off
xlabel("MOG eccentricity");
ylabel("Slope of resulting distribution [km/s/rad]");
legend("data", "fit")


a_coeffs = zeros(size(as));
b_coeffs = zeros(size(as));
c_coeffs = zeros(size(as));
% d_coeffs = zeros(size(as));
RMSEs = zeros(size(as));

for i = 1:length(as)
    cs = coeffvalues(fits.("fit"+i));
    a_coeffs(i) = cs(1);
    b_coeffs(i) = cs(2);
    c_coeffs(i) = cs(3);
    % d_coeffs(i) = cs(4);
    RMSEs(i) = gofs.("gof" + i).rmse;
end

figure(27949574)
plot(as/AU, a_coeffs)
xlabel("aMOG [AU]");
ylabel("a coefficient");

figure(27949575)
plot(as/AU, b_coeffs)
xlabel("aMOG [AU]");
ylabel("b coefficient");

figure(27949576)
plot(as/AU, RMSEs)
xlabel("aMOG [AU]");
ylabel("RMSE");

figure(27949577)
plot(as/AU, RMSEs./(mean(coeffs, 2))')
xlabel("aMOG [AU]");
ylabel("RMSE [%]");

%% Construct 2D Fit

% b = mean(b_coeffs);
[fita, gofa] = fit(as'/AU, a_coeffs', "a*x^(-0.5)");
[fitb, gofb] = fit(as'/AU, b_coeffs', "a*x^(-0.5)");
[fitc, gofc] = fit(as'/AU, c_coeffs', "a*x^(-0.5)");
% [fitd, gofd] = fit(as'/AU, d_coeffs', "a*x^(-0.5)");

%% Evaluating numeric fit accuracy
% See MOGPhasingNumericFit function

% coeffs_numeric = MOGPhasingNumericFitPositive(Y*AU, X);
coeffs_numeric = MOGPhasingNumericFitNegative(Y*AU, X);

errors = coeffs_numeric - coeffs;
errors_percent = (errors)./coeffs * 100;

figure(28392830)
surf(X,Y,errors)
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]")
zlabel("Error from numeric approximation [km/s/rad]")

figure(28392831)
surf(X,Y,errors_percent)
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]")
zlabel("Error from numeric approximation [%]")

es2 = 1e-2:1e-3:5e-1;
as2 = AU*(0.7:0.01:1.5);
[X2,Y2] = meshgrid(es2, as2/AU);
% coeffs_numeric2 = MOGPhasingNumericFitPositive(Y2*AU, X2);
coeffs_numeric2 = MOGPhasingNumericFitNegative(Y2*AU, X2);

figure(28392829)
surf(X2,Y2,coeffs_numeric2)
hold on
surf(X,Y,coeffs)
hold off
xlabel("MOG eccentricity");
ylabel("MOG semi-major axis [AU]")
zlabel("Slope of numerically-approximated distribution [km/s/rad]");

%% TOTAL RMSE

RMSE = rmse(coeffs_numeric, coeffs, "all")
RMSE_percent = RMSE/mean(coeffs,"all") * 100
MAE = mae(coeffs_numeric, coeffs, "all")
MAE_percent = MAE/mean(coeffs,"all") * 100


%% Save Data

% savefile = "../Plotting/PlottingData/MOG Phasing/MOGSurfaceData1Neg.mat";
% save(savefile)

%% Functions

function coeff = linearfit(a,e,direction)
    if direction == 1
        dPs = linspace(0.0, 2*pi/20, 5);
        % dPs = linspace(0.0, -1, 20); % Testing GOF
    elseif direction == -1
        dPs = linspace(0.0, -2*pi/20, 5);
    end
    DVs = Inf(size(dPs));
    ToPs = DVs;
    muSu = 1.327124400419393e+11;
    
    for i = 1:length(dPs)
        dP = dPs(i);
        [DV, ToP] = computeOptimalPhasingDVMOG(a, e, dP, muSu);
        DVs(i) = DV;
        ToPs(i) = ToP;
    end
    options = fitoptions( 'Method', 'NonlinearLeastSquares',...
    'Lower',[-100], 'Upper',[100], 'StartPoint', 5);
    [fit1, gof1] = fit(dPs(1:end)', DVs(1:end)', 'a*x', options);
    coeff = coeffvalues(fit1);
end