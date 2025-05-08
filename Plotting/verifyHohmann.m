% Compare Approaches

AU = 149600000; % km
muSu = 1.327124400419393e+11; 
NA = 5;
NE = 6;

aMOGs = linspace(0.7, 1.5, NA) * AU; % km
eMOGs = linspace(0.1, 0.6, NE);

Slot_DVs = NaN([NA, NE]);
Hohmann_DVs = NaN([NA, NE]);

%% Run

for i = 1:NA
    aMOG = aMOGs(i);
    for j = 1:NE
        eMOG = eMOGs(j);

        [MOGSat_DV, ~, ~] = findOptimalTFC2MOGSlot(aMOG,eMOG,muSu);
        [DV, ~, ~, ~, ~] = computeHohmannTransferCirc2Ell(aMOG,aMOG,eMOG,muSu);

        Slot_DVs(i,j) = MOGSat_DV;
        Hohmann_DVs(i,j) = DV;

        disp("Progress: " + num2str((i-1)*NE + j) + "/" + num2str(NE*NA) + " analyses done")
    end
end

%% Analysis

difference = Slot_DVs - Hohmann_DVs;
difference>=0