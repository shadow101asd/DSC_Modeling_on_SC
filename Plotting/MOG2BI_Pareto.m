
AU = 149600000;
aMOG = AU;
day = 24*3600;
% Ne = 9;
% eMOGs = linspace(1/(Ne+1), 1-1/(Ne+1), Ne);
eMOGs = 0.05:0.05:0.95;
Ne = length(eMOGs);
muSu = 1.327124400419393e+11;

res = 1000;

%% Run

nT = 100;

TOFs = [];
DVs = NaN([Ne, nT]);
Ki = [aMOG, 0, 0, 0, 0, 0];

for i = 1:Ne
    eMOG = eMOGs(i);
    Kf = [aMOG, eMOG, 0, 0, pi, 0];
    
    TOF_hohmann = pi*sqrt((aMOG*(1+1/2 *eMOG))^3 /muSu);
    TOFs = [TOFs; linspace(0, TOF_hohmann, nT)];

    for j = 1:nT
        TOF = TOFs(i,j);
        if j < 0.2*nT
            DVs(i,j) = findMinDVTransferForFixedTOF(Ki, Kf, TOF, 2.0*res, muSu); % Increase resolution near edge for numerical stability
        else
            DVs(i,j) = findMinDVTransferForFixedTOF(Ki, Kf, TOF, res, muSu);
        end

        disp("Progress: " + num2str((i-1)*nT+j) + "/" + num2str(Ne*nT) + " operations")
    end
end

%% Analysis

% Local Minimum Solution

TF = islocalmin(DVs, 2);
LocalMin_DVs_temp =  DVs(TF)';
LocalMin_TOFs_temp = TOFs(TF)';

% Remove outliers
[~,maxidx] = max(LocalMin_DVs_temp);
LocalMin_DVs_temp = LocalMin_DVs_temp(maxidx:end);
LocalMin_TOFs_temp = LocalMin_TOFs_temp(maxidx:end);
LM_start = length(eMOGs) - length(LocalMin_DVs_temp) + 1; % Index of first local min.

% Reformat
LocalMin_DVs = NaN(size(eMOGs));
LocalMin_DVs(LM_start:end) = flip(LocalMin_DVs_temp);
LocalMin_TOFs = NaN(size(eMOGs));
LocalMin_TOFs(LM_start:end) = flip(LocalMin_TOFs_temp);

clear LocalMin_TOFs_temp LocalMin_DVs_temp

% Hohmann Calcs

Hohmann_DVs = zeros([1,Ne]);
Hohmann_TOFs = zeros([1,Ne]);
for i = 1:Ne
    [DV, ToF] = computeHohmannTransferCirc2Ell(aMOG,aMOG,eMOGs(i),muSu);
    Hohmann_DVs(i) = DV;
    Hohmann_TOFs(i) = ToF;
end

% 1B Insertion Calcs
OneB_DVs = sqrt(muSu/aMOG) * sqrt(2*(1-sqrt(1-eMOGs.^2)));

%% Plotting

figure(293239)
plot(TOFs'/day, DVs')
xlabel("Time of flight [days]")
ylabel("2-burn Transfer DV [km/s]")
legend(string(eMOGs))

% AAS MOG Paper Plotting

figure(129812)
plot(eMOGs',OneB_DVs)
hold on
plot(eMOGs', LocalMin_DVs)
plot(eMOGs',Hohmann_DVs)
hold off
xlabel("eMOG")
ylabel("Transfer DV [km/s]")
legend("1-burn Transfer", "Local Min. Transfer", "Hohmann Transfer")

figure(129813)
plot(eMOGs',0*Hohmann_TOFs/day) % One-burn transfer has zero flight time
hold on
plot(eMOGs', LocalMin_TOFs/day)
plot(eMOGs',Hohmann_TOFs/day)
hold off
xlabel("eMOG")
ylabel("Transfer Time of Flight (TOF) [days]")
legend("1-burn Transfer", "Local Min. Transfer", "Hohmann Transfer")