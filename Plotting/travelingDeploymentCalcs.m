function [Shuttle_DVs, Sat_DVs] = travelingDeploymentCalcs(aMOG, eMOG, mu, deployment_times, direction)
% Docstring

day = 3600*24; % a day in seconds

K0 = [aMOG; eMOG; 0; 0; -pi/2; pi];
N = length(deployment_times);
dPs = linspace(0, direction*2*pi, N);

Shuttle_DVs = zeros(size(deployment_times));
Sat_DVs = Shuttle_DVs;

V1s = zeros([3, N]);
V2s = V1s;

for i = 1:(N-1)
    t1 = deployment_times(i);
    t2 = deployment_times(i+1);

    K1 = shiftK_inMOG(K0, dPs(i), mu);
    K2 = shiftK_inMOG(K0, dPs(i+1), mu);

    X1 = propagateFromKeplerians(K1, mu, [0, t1]);
    X2 = propagateFromKeplerians(K2, mu, [0, t2]); 

    X1 = X1(:,2);
    X2 = X2(:,2);

    R1 = X1(1:3)';
    R2 = X2(1:3)';
    
    % Lambert
    tf = (t2-t1)/day; % days

    [theta1, ~] = cart2pol(X1(1), X1(2));
    [theta2, ~] = cart2pol(X2(1), X2(2));
    dtheta = mod(theta2-theta1, 2*pi);
    if dtheta > pi
        tf = -tf; % Take the long path in the lambert solver
    end

    m = 0; % No loops

    [V1, V2, ~, exitflag] = lambert(R1, R2, tf, m, mu);
    V1s(:,i) = V1';
    V2s(:,i) = V2';

    % Compute Shuttle DVs
    if i == 1 %% Shuttle DV is from MOG velocity at the first step
        SDV = abs(norm(X1(4:6)'-V1));
    else
        SDV = abs(norm(V2s(:,i-1)'-V1)); % Compare Incoming and Outgoing DV at each step
    end
    Shuttle_DVs(i) = SDV;

    % Compute SatDVs
    % Sat could be dropped off from arrival, or departure trajectory, so
    % we'll consider both cases

    if i == 1 
        SATDV = 0; % We're assuming that the Shuttle starts in the MOG, so we get the first dropoff for "free"
    else
        arrivalV = V2s(:,i-1);
        departureV = V1s(:,i);

        MOGV = X1(4:6);
        
        DVA = abs(norm(arrivalV-MOGV));
        DVB = abs(norm(departureV-MOGV));

        SATDV = min(DVA,DVB);
    end
    
    Sat_DVs(i) = SATDV;
    
end

% Compute the last SatDV:

t1 = deployment_times(N);
K1 = shiftK_inMOG(K0, dPs(N), mu);
X1 = propagateFromKeplerians(K1, mu, [0, t1]);
X1 = X1(:,2);

arrivalV = V2s(:,N-1);
departureV = V1s(:,N);

MOGV = X1(4:6);

DVA = abs(norm(arrivalV-MOGV));
DVB = abs(norm(departureV-MOGV));

Sat_DVs(N) = min(DVA,DVB);

end