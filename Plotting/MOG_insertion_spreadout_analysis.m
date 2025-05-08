%% Parabola comp.
AU = 149600000; % AU in km

NSats = 20;
eMOG = 0.3;
aMOG = 1.0*AU;
muSu = 1.327124400419393e+11;
Sat_Isp = 300; % s
Total_Payload_M = 200*10^3; % Max Starship Block III Payload

[m_persat, maxTOP] = computeSatMForPhasingDeployment(aMOG, eMOG, muSu, NSats, Sat_Isp, Total_Payload_M)
