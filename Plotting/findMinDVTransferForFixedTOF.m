function [best_DV, best_X1, best_X2, best_DV1, best_DV2] = findMinDVTransferForFixedTOF(Ki, Kf, TOF, res, mu)
    day = 3600*24; % a day in seconds
    
    % Orbital periods
    T_i = 2*pi*sqrt(Ki(1)^3/mu);
    T_f = 2*pi*sqrt(Kf(1)^3/mu);
    etR_i = linspace(0, T_i, res);
    etR_f = linspace(0, T_f, res);
    Ni = length(etR_i);
    Nf = length(etR_f);
    
    % Generate Orbits
    Xi = propagateFromKeplerians(Ki, mu, etR_i);
    Xf = propagateFromKeplerians(Kf, mu, etR_f);

    % Initialize variables
    best_DV = inf;
    best_X1 = [];
    best_X2 = [];
    best_DV1 = inf;
    best_DV2 = inf;

    for i = 1:Ni
        for f = 1:Nf
            Ri = Xi(1:3, i)';
            Rf = Xf(1:3, f)';

            % Check to see if TOF needs to be negated
            [thetai, ~] = cart2pol(Ri(1), Ri(2));
            [thetaf, ~] = cart2pol(Rf(1), Rf(2));
            dtheta = angdiff(thetai, thetaf); % in [-pi, pi]
            
            if dtheta >= 0 % Enforcing CCW transfers only
                [V1, V2, ~, exitflag] = lambert(Ri, Rf, TOF/day, 0, mu);
            else % Enforcing CCW transfers only
                [V1, V2, ~, exitflag] = lambert(Ri, Rf, -TOF/day, 0, mu);
            end
            
            if exitflag
                Vi = Xi(4:6, i)';
                Vf = Xf(4:6, f)';
                DV1 = norm(V1 - Vi);
                DV2 = norm(Vf - V2);
                DV = abs(DV1) + abs(DV2);

                if DV < best_DV
                    best_DV = DV;
                    best_DV1 = DV1;
                    best_DV2 = DV2;
                    best_X1 = Xi(:, i);
                    best_X2 = Xf(:, f);
                end
            end
        end
    end


end