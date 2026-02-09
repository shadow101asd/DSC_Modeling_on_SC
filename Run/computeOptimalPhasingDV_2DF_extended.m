function [best_DV, ToP, DV1, DV2, details] = computeOptimalPhasingDV_2DF_extended(W, M, eF, dPhase)
    mu = 1.327124400419393e+11; % Sun gravitational parameter
    AU = 149600000; % 1 AU in km
    aF = AU;

    % Reduce W,M to fundamental pair
    G = gcd(W,M);
    W = W/G;
    M = M/G;

    % Meta-parameters
    day = 86400; % 1 day in seconds
    max_periods = 2.0; % in R+ | TODO: adapt to represent full range of possible starting positions
    T = 2*pi*sqrt(aF^3/mu); % satellites' orbital period
    min_tf_considered = 0.2*T/day; % in days
    max_tf_considered = 2.0*T/day; % in days
    m = 0;
    
    N = 150; % desired # of sample points
    dt = day*ceil(T*max_periods/day/N); % s THIS CAN BE LARGER NOW WITH THE EXTENSION

    % Setup
    et0 = 0; % Arbitrary sim. start time
    et1 = et0 + max_periods*T;
    etR = et0:dt:et1;

    K1 = [aF; eF; 0; 0; -pi/2; pi];
    K2 = shiftK_in_2DF(K1, W, M, dPhase, mu);

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

            theta1 = atan2(X1(2),X1(1)); % The exact same as cart2pol since we don't care about r
            theta2 = atan2(X2(2),X2(1));

            dtheta = mod(theta2-theta1, 2*pi);
            if dtheta > pi
                tf = -tf; % Take the long path in the lambert solver
            end
            
            R1 = X1(1:3)';
            R2 = X2(1:3)';
            
            try
                [V1, V2, ~, exitflag] = lambert_mex(R1, R2, tf, m, mu);
            catch ME
                % use uncompiled version if the compiled version fails (or
                % doesn't exist?)
                [V1, V2, ~, exitflag] = lambert(R1, R2, tf, m, mu); 
            end

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

    % EXTENSION: After Global best is found, run fmincon to find local
    % optimum
    
    factor = 1e6;
    x0 = [details.t1, details.t2]/factor;
    A = [1 -1; 0 0];
    b = [0; 0];
    lb = [0; 0];
    ub = [et1-min_tf_considered*day ; et1]/factor;
    options = optimoptions("fmincon", "Display", "none", OptimalityTolerance=1e-4);

    [x_opt,fval] = fmincon(@(X) lambertWrapper(X, K1, K2, mu, factor), x0, A, b, [], [], lb, ub, [], options);

    % Replace Outputs with new best

    best_t1 = x_opt(1)*factor;
    best_t2 = x_opt(2)*factor;
    ToP = best_t2-best_t1;
    [best_DV, DV1, DV2, details_L] = lambertWrapper(x_opt, K1, K2, mu, factor);
    details.t1 = best_t1;
    details.t2 = best_t2;
    details.toP_days = (best_t2-best_t1)/day;
    details.V1 = details_L.V1;
    details.V2 = details_L.V2;
    details.R1 = details_L.R1;
    details.R2 = details_L.R2;
    

    function [DV, DV1, DV2, details] = lambertWrapper(X, K1, K2, mu, factor)
        t1 = X(1)*factor;
        t2 = X(2)*factor;
        tf_w = (t2 - t1)/day;

        X1s_w = propagateFromKeplerians(K1, mu, [0, t1]);
        X2s_w = propagateFromKeplerians(K2, mu, [0, t2]);

        X1_w = X1s_w(:,2);
        X2_w = X2s_w(:,2);

        [theta1_w, ~] = cart2pol(X1_w(1), X1_w(2));
        [theta2_w, ~] = cart2pol(X2_w(1), X2_w(2));
        dtheta_w = mod(theta2_w-theta1_w, 2*pi);
        if dtheta_w > pi
            tf_w = -tf_w; % Take the long path in the lambert solver
        end
            
        R1_w = X1_w(1:3)';
        R2_w = X2_w(1:3)';
        [V1_w, V2_w, ~, exitflag_w] = lambert(R1_w, R2_w, tf_w, 0, mu);

        if exitflag_w == 1
           DV1 = abs(norm(V1_w-X1_w(4:6)'));
           DV2 = abs(norm(V2_w-X2_w(4:6)'));
    
           DV = DV1 + DV2;
           details.V1 = V1_w;
           details.V2 = V2_w;
           details.R1 = R1_w;
           details.R2 = R2_w;
        end
        
    end
end