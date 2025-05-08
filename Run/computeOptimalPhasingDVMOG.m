function [best_DV, ToP, DV1, DV2, details] = computeOptimalPhasingDVMOG(aMOG, eMOG, dPhase, mu)
    % Meta-parameters
    day = 86400; % 1 day in seconds
    temporal_res = 100; % per TMOG TO ADJUST FOR DECENT CONVERGENCE
    max_periods = 2.0; % in R+
    TMOG = 2*pi*sqrt(aMOG^3/mu); % MOG orbital period
    dt = TMOG/temporal_res; % s
    min_tf_considered = 0.2*TMOG/day; % in days
    max_tf_considered = 0.9*TMOG/day; % in days
    m = 0;

    % Setup
    et0 = 0; % Arbitrary sim. start time
    et1 = et0 + max_periods*TMOG;
    etR = et0:dt:et1;

    K1 = [aMOG; eMOG; 0; 0; -pi/2; pi];
    K2 = shiftK_inMOG(K1, dPhase, mu);

    X1s = propagateFromKeplerians(K1, mu, etR);
    X2s = propagateFromKeplerians(K2, mu, etR);
    
    nT = length(etR);
    best_DV = Inf;
    details = struct;
    details.t1 = -1;
    details.t2 = -1;
    for i = 1:nT
        for j = i:nT
            tf = (j-i)*dt/day; % tof in days

            if tf < min_tf_considered || tf > max_tf_considered
                continue % No need to consider this pairing
            end

            X1 = X1s(:,i);
            X2 = X2s(:,j);

            [theta1, ~] = cart2pol(X1(1), X1(2));
            [theta2, ~] = cart2pol(X2(1), X2(2));
            dtheta = mod(theta2-theta1, 2*pi);
            if dtheta > pi
                tf = -tf; % Take the long path in the lambert solver
            end
            
            R1 = X1(1:3)';
            R2 = X2(1:3)';
            [V1, V2, ~, exitflag] = lambert(R1, R2, tf, m, mu);

            if exitflag == 1
                DV1_ij = abs(norm(V1-X1(4:6)'));
                DV2_ij = abs(norm(V2-X2(4:6)'));
    
                DV = DV1_ij + DV2_ij;
    
                if DV < best_DV
                    best_DV = DV;
                    DV1 = DV1_ij;
                    DV2 = DV2_ij;
                    ToP = abs(tf)*day;
                    details.t1 = etR(i);
                    details.t2 = etR(j);
                    details.toP_days = abs(tf);
                end
            end

        end
    end

    
end