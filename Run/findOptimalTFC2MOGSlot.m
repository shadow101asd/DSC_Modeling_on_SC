function [MOGSat_DV, LToF, lambert_details] = findOptimalTFC2MOGSlot(aMOG,eMOG,mu)
% findOptimalTFC2MOGSlot Summary goes here: Could be further vectorized?
% THIS IS OBSOLETE! USE computeHohmannTransferCirc2Ell instead
tic
% Optimization Meta-Parameters
NMOG_candidates = 100;
day = 86400; % 1 day in seconds
dt = day*3; % s
max_periods = 0.9; % in R+
m_max = 0;

% Setup
et0 = 0; % Arbitrary sim. start time
TMOG = 2*pi*sqrt(aMOG^3/mu); % MOG orbital period
et1 = et0 + max_periods*TMOG;
etR = et0:dt:et1; % Row vector of time of flights
et2 = et0 + max_periods*TMOG;
etR_1rev = et0:dt:et2; % Row vector of times for candidate slot positions

% Generate slot positions over one full period
XP = Keplerian2Cartesian(aMOG, 0, 0, 0, 0, pi/2, mu);
XSlots = generateMOGNearPlanet(XP, mu, etR_1rev, eMOG, 1, NMOG_candidates,false);

% Initial conditions
R1 = [0, aMOG, 0];
V1_initial = [-sqrt(mu/aMOG), 0, 0];

% Data arrays
minDVs = Inf(size(etR));
minIdxs = -ones(2, size(etR, 2));
ms = -NaN(size(etR));

% Arrival time loop
for tf_idx = floor(TMOG/(4*dt)):length(etR)
    t = etR(tf_idx);
    tf = (t-et0)/day; % Time of flight, in days. Positive for the short route around the Sun!
    % m_max = floor((t-et0)/TMOG*1.2); % Number of revolutions of the MOG around the Sun before arrival
    
    % Slot candidates loop
    for i = 1:NMOG_candidates
        % Positions during one MOG orbit
        for j = 1:length(etR_1rev)
            % j = max([mod(tf_idx, length(etR_1rev)),1]); % Testing to
            % reobtain previous, fixed phased results

            % Arrival conditions
            R2 = XSlots(1:3,j,i)';
            V2_final = XSlots(4:6,j,i)';
            
            % M loop
            for m = -m_max:m_max
                % Run lambert solver
                if R2(1) <= 0
                    [V1, V2, ~, exitflag] = lambert(R1, R2, tf, m, mu);
                else
                    % Take the "long way around" if CCW is the long way around
                    [V1, V2, ~, exitflag] = lambert(R1, R2, -tf, m, mu);
                end
        
                if exitflag % Solver exited correctly
                    DV1 = norm(V1-V1_initial);
                    DV2 = norm(V2_final-V2);
                    DV = DV1 + DV2;
        
                    if DV < minDVs(tf_idx) % Found a better solution for that timestep
                        minDVs(tf_idx) = DV;
                        minIdxs(1,tf_idx) = i;
                        minIdxs(2,tf_idx) = j;
                        ms(tf_idx) = m;
                    end
                end
            end
        end
    end
end

toc

[MOGSat_DV, min_idx] = min(minDVs); 
LToF = etR(min_idx)-et0; % Time of flight, in s

lambert_details = []; % Can add other data later on.

end