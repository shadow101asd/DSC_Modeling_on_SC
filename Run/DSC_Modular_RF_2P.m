function [] = DSC_Modular_RF_2P(Nl, Nf, run_idx, gaoptions, sat_config_name, shuttle_config_name, plotting_bool, warm_start_bool, Planet_string, fmincon_bool)
    
    Nf_original = Nf; % make a copy for file saving
    disp("Nl = " + Nl) % Print number of launches, for tracking
    if Nf <= Nl
        disp("Nf = " + Nf) % Print number of features/blocks, for tracking
    else
        disp("Nf > Nl. Fixing Nf = Nl = " + Nl) % Print number of features/blocks, for tracking
    end

    % Constants

    AU = 1.496e8; % 1 AU in km
    
    % Load ephemerides and other pregenerated data.
    
    if nargin <= 8
        Planet_string = 'Mars';
        disp("No Planet string passed in, defaulting to Mars.")
    else
        disp("Planet: " + Planet_string)
    end
    

    load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "etR");
    PX_string = "X" + Planet_string(1:2);
    XP2 = load("../Inputs/run"+run_idx+".mat", PX_string).(PX_string);

    % Load and process satellite and terminal specs
    load("../Inputs/SE/"+sat_config_name+".mat", "Sat_specs", "Earth_specs", "Comms_specs");
    P2_specs = load("../Inputs/SE/"+sat_config_name+".mat", Planet_string + "_specs").(Planet_string + "_specs");
    P2_specs.pstring = Planet_string;

    % Load Starship specs
    load("../Inputs/SE/"+shuttle_config_name+".mat", "Shuttle_specs");
    
    % Compute bounds on semi-major axes of features
    as = [mean(Cartesian2Keplerian(XEa, muSu)), mean(Cartesian2Keplerian(XP2, muSu))];
    min_a = min(as)/AU; 
    max_a = max(as)/AU; 

    % Can't precompute link budget matrix due to variable number of
    % satellites

    % Set up optimization
    lb_ef     = [1 0  1 min_a 0   1e-6	     1e-6      0 1];
    ub_ef     = [6 Nl 2 max_a 0.7 2*pi-1e-6  2*pi-1e-6 1 1];
    % ub_ef     = [8 Nl 2 max_a 1 2*pi 1 5];
    intcon_ef = [1 2 3 9];
    
    nvars_pf = 9;
    nvars = nvars_pf*Nf;

    lb = repmat(lb_ef, [1, Nf]);
    ub = repmat(ub_ef, [1, Nf]);

    if Nf > Nl
        ub(2+Nl*nvars_pf:nvars_pf:end) = 0; % fix later features to 0 launches to reduce search space
    end

    lb(2) = 1; % Enforce at least one launch in the entire constellation

    intcon = zeros(1,4*Nf);
    for i = 1:Nf
        intcon(1+4*(i-1):4*i) = intcon_ef + (i-1)*nvars_pf;
    end
    
    A = zeros([Nf, nvars]);
    b = zeros([Nf, 1]);
    
    % Remove symmetries to shrink design space by a factor of 1/Nf factorial
    for i = 1:(Nf-1)
        A(i, 1 + (i-1)*nvars_pf) = 1;
        A(i, 1 + i*nvars_pf) = -1;
    end
    
    % Enforce total # of launches less than Nl input
    A(end, 2:nvars_pf:end) = 1;
    b(end) = Nl;


    nonlcon = @(x) deal(A*x' - b, []);
    % nonlcon = []; Testing
  
    % A = [];
    % b = []; % Testing
    

    % Tweak gaoptions according to inputs and Nf, Nl
    
    if plotting_bool % Tie in appropriate custom plotting function
        plotFcnWrapper = @(options, state, flag) ...
                        plotBestConstellation(options, state, flag, ...
                          XEa, XP2, muSu, Shuttle_specs, ...
                          Comms_specs, Sat_specs, Earth_specs, ...
                          P2_specs, Nf);
        gaoptions = optimoptions(gaoptions, 'PlotFcn', plotFcnWrapper);
    end


    % Modify gaoptions Population Count to scale affinely with Nf
    gaoptions.PopulationSize = floor(gaoptions.PopulationSize * (1 + (Nf-1)/2)); 

    % Include warmstart if requested
    if warm_start_bool
        N_warmstart = floor(0.1*gaoptions.PopulationSize);
        Initial_Pop = generate_WarmStart_Individuals(N_warmstart, run_idx, Nf, Nl, [ub_ef(1),lb_ef(2:end)], Nf);

        gaoptions.InitialPopulationMatrix = Initial_Pop;
        disp(int2str(size(Initial_Pop, 1))+"/"+int2str(N_warmstart)+ " individuals found for the warm start.") % for monitoring
    end


    % Run the Genetic Algorithm
    func = @(X) wrapperFunc_RF_2P(X, XEa, XP2, etR, muSu, Shuttle_specs, Comms_specs, Sat_specs, Earth_specs, P2_specs, Nf);
    
    [X_opt, fval, EXIT_FLAG, output] = ga(func, nvars,[],[],[],[],lb,ub,nonlcon,intcon,gaoptions);
    
    disp("GA Output: " + num2str(-fval)) % Show solution in logs

    % Hybrid function
    if nargin <= 9
        fmincon_bool = true; % Default to true
    end
    if fmincon_bool
        disp("Finishing with hybrid_fcn: fmincon");
        
        % Modify lb and ub to fix integer variables
        lb_intsRfixed = lb;
        lb_intsRfixed(1:nvars_pf:end) = X_opt(1:nvars_pf:end);
        lb_intsRfixed(2:nvars_pf:end) = X_opt(2:nvars_pf:end);
        lb_intsRfixed(3:nvars_pf:end) = X_opt(3:nvars_pf:end);
        lb_intsRfixed(9:nvars_pf:end) = X_opt(9:nvars_pf:end);

        ub_intsRfixed = ub;
        ub_intsRfixed(1:nvars_pf:end) = X_opt(1:nvars_pf:end);
        ub_intsRfixed(2:nvars_pf:end) = X_opt(2:nvars_pf:end);
        ub_intsRfixed(3:nvars_pf:end) = X_opt(3:nvars_pf:end);
        ub_intsRfixed(9:nvars_pf:end) = X_opt(9:nvars_pf:end);

        % Run fmincon
        fmincon_opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'UseParallel', true);
        
        if plotting_bool % Tie in appropriate plotting function for fmincon
            fmincon_opts = optimoptions(fmincon_opts, 'PlotFcn', 'optimplotfval');
        end

        [X_opt, fval_f] = fmincon(func, X_opt, [], [], [], [], lb_intsRfixed, ub_intsRfixed, [], fmincon_opts);

        disp("Fmincon Output: " + num2str(-fval_f)) % Show fval in logs
    end

    % Save data

    [Out, XSats, NSats] = wrapperFunc_RF_2P(X_opt, XEa, XP2, etR, muSu, Shuttle_specs, Comms_specs, Sat_specs, Earth_specs, P2_specs, Nf);
    
    B = -Out;

    try 
        XSats_i = XSats(:,1,:); % Initial Satellite Positions and Velocities
    catch ME
        XSats_i = [];
        warning("Optimal solution includes no satellites.")
    end
    
    filename = "../Data/run"+run_idx+"/Nf"+int2str(Nf_original)+"Nl"+int2str(Nl)+".mat";

    if isfile(filename)
        old_B = load(filename, "B").B;
        if old_B < B
            save(filename, "B", "X_opt", "NSats", "XSats_i", "output", "lb", "ub");
        end
    else
        save(filename, "B", "X_opt", "NSats", "XSats_i", "output", "lb", "ub");
    end
    
    % Final message.
    disp("Final solution: " + num2str(X_opt))
    disp("evaluating to a bandwidth of " + num2str(B) + " Mbps.")

end

function [Out, XSats, NSats] = wrapperFunc_RF_2P(X, X1, X2, etR, mu, Shuttle_specs, Comms_specs, Sat_specs, P1_specs, P2_specs, Nf)
    
    % Reshape X for convenience
    X = reshape(X, 9, Nf)';

    % Compute NSats for each feature

    NSats = zeros([Nf, 1]);
    D_antennas = zeros([Nf, 1]);
    for i = 1:Nf
        [a, b] = getNSats(X(i,:), {X1, X2}, mu, Shuttle_specs, Sat_specs);
        NSats(i) = a;
        D_antennas(i) = b;
    end

    % Compute Link Budget Matrix
    R_1km = buildAdjacency_LinkBudgetMatrix(NSats, Comms_specs, Sat_specs, P1_specs, P2_specs, D_antennas); % Units: bps*km^2

    % Compute XSats
    
    XSats = [];
    for i = 1:Nf
        XSats_i = getXSats(X(i,:), {X1, X2}, mu, etR, NSats(i));
        if isempty(XSats)
            XSats = XSats_i;
        else
            XSats = cat(3,XSats, XSats_i);
        end
    end

    % Evaluate Bandwidth (pass in max datarate due to hardware limitations)
    Out = bestLinkBudget_bandwidth(X1,X2,XSats,R_1km,Comms_specs.maxDR_Mbps,P2_specs.pstring);
end

function [Block_ID, N_launches, Planet_ID, a, e, w, f0, frac, Npl] = unpackVars(X)
    AU = 149600000; % 1AU in km

    Block_ID = round(X(1));
    N_launches = round(X(2));
    Planet_ID = round(X(3));
    a = X(4) * AU; % convert to km
    e = X(5);
    w = X(6);
    f0 = X(7);
    frac = X(8);
    Npl = round(X(9));

    if Block_ID == 7 || Block_ID == 8
        N_launches = max(N_launches, 1);
    end
end

function [NSats, D_antennas] = getNSats(X, XPs, mu, Shuttle_specs, Sat_specs)
    maxO = 3;

    [Block_ID, N_launches, Planet_ID, a, e, ~, ~, frac, Npl] = unpackVars(X);

    if N_launches == 0 % We can end early
        NSats = 0;
        D_antennas = 0;
        return
    end
    
    % Define function for finding mpersat as a function of Nsats only, for
    % each case

    switch(Block_ID)
        case 1 % Regular MOG
            frac_pl = 1/N_launches;
            mps_func = @(N) getMOG_SatDryMasses(a,e,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, frac_pl);

        case 2 % MOG @ Planet

            % Get Corresponding Planet Info
            XP = XPs{Planet_ID};
            aP = mean(Cartesian2Keplerian(XP,mu));

            frac_pl = 1/N_launches;
            mps_func = @(N) getMOG_SatDryMasses(aP,e,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, frac_pl);

        case 3 % MOG Including Planet

            % Get Corresponding Planet Info
            XP = XPs{Planet_ID};
            [aP, eP] = Cartesian2Keplerian(XP,mu);
            aP = mean(aP);
            eP = mean(eP);

            frac_pl = 1/N_launches;
            mps_func = @(N) getMOG_SatDryMasses(aP,eP,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, frac_pl);

        case 4 % Circular Ring

            frac_pl = frac / N_launches;
            mps_func = @(N) getSatDryMasses_CircularLoop(a,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, maxO, frac_pl);

        case 5 % Elliptical Ring

            frac_pl = frac / N_launches;
            mps_func = @(N) getSatDryMasses_EllipticalLoop(a,e,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, maxO, frac_pl);
            
        case 6 % Ring Along a Planet's Orbit

            % Get Corresponding Planet Info
            XP = XPs{Planet_ID};
            [aP, eP] = Cartesian2Keplerian(XP,mu);
            aP = mean(aP);
            eP = mean(eP);

            frac_pl = frac / N_launches;
            mps_func = @(N) getSatDryMasses_EllipticalLoop(aP,eP,mu,N,Shuttle_specs.Isp,Sat_specs.Isp,Shuttle_specs.wetMass,Shuttle_specs.dryMass,Shuttle_specs.maxPayload, maxO, frac_pl);

        case 7 % Satellite @ Lagrange Point

            % TODO: make functions for these
            error("Not yet implemented.")

        case 8 % Lone Satellite
            
            % TODO: make functions for these
            f0 = frac*2*pi; % Rescale to correct domain
            error("Not yet implemented.")
             
        otherwise
            error("Incorrect Building Block ID")
    end
    
    % Find Nsats per launch
    [NSats_pl, D_antennas] = compute_Nsats_per_Shuttle_RF(Shuttle_specs, Sat_specs, Npl, mps_func);

    % Find total Nsats in the block
    NSats = NSats_pl * N_launches;

    if isnan(NSats) % Clean up output
        NSats = 0;
        D_antennas = 0;
    end

end



function XSats = getXSats(X, XPs, mu, etR, Nsats)
   
    % Generate Cartesian Coordinates of the satellites described by this
    % row of the genome

    [Block_ID, ~, Planet_ID, a, e, w, f0, frac, Npl] = unpackVars(X);

    switch(Block_ID)
        case 1 % Regular MOG

            Ki = [a; 0.25; 0.0; 0.0; w; 0.0]; % Second element (eccentricity) isn't actually used here
            XSats = NSATSpropagateFromKepleriansSHELLS(Ki,mu,etR,Nsats,e,1);

        case 2 % MOG @ Planet

            XSats = generateMOGNearPlanet(XPs{Planet_ID}, mu, etR, e, 1, Nsats, true);
            
        case 3 % MOG Including Planet

            XSats = generateMOGIncludingPlanet(XPs{Planet_ID}, mu, etR, Nsats);

        case 4 % Circular Ring

            XSats = NSATSpropagateFromKepleriansCIRC(a,f0,mu,etR,Nsats,frac);

        case 5 % Elliptical Ring TODO
            
            Xcenter = Keplerian2Cartesian(a,e,0,0,w,f0,mu);
            XSats = generateXsAlongOrbit(Xcenter, mu, etR, Nsats, frac);
            
        case 6 % Ring Along a Planet's Orbit

            XSats = generateXsAlongOrbit(XPs{Planet_ID}, mu, etR, Nsats, frac);

        case 7 % Satellite @ Lagrange Point TODO

            L4_L5 = round(frac);

        case 8 % Lone Satellite TODO

            
             
        otherwise
            error("Incorrect Building Block ID")
    end

end

function state = plotBestConstellation(options, state, flag, X1, X2, mu, Shuttle_specs, Comms_specs, Sat_specs, P1_specs, P2_specs, Nf)
    % Called by ga() at each generation

    % h = 420; % Figure number

    % Reuse GA figure
    h = findall(0, 'Type', 'figure', 'Name', 'Genetic Algorithm');
    if isempty(h)
        h = figure('Name', 'Genetic Algorithm');
    else
        figure(h); clf;
    end

    % Get best individual

    [B,best_idx] = min(state.Score);
    bestGenome = state.Population(best_idx,:);  % the best of current gen
            
    % Compute XSats and associated metrics for plotting
    
    dummy_etR = 1;

    [~, XSats, NSats] = wrapperFunc_RF_2P(bestGenome, X1(:,1), X2(:,1), dummy_etR, mu, Shuttle_specs, Comms_specs, Sat_specs, P1_specs, P2_specs, Nf);
    XSats_flat = reshape(XSats, [6,sum(NSats)]);
    AU = 149600000; % 1 AU in km

    % Plot or update your custom plot
    figure(h);  % Bring focus to figure
    
    scatter(0, 0, "yellow", 'filled', 'hexagram', MarkerEdgeColor='black');
    hold on
    scatter(X1(1,1)/AU, X1(2,1)/AU, "blue", 'filled');

    % Break up plotting cases depending on the specific planet
    if  strcmp(string(P2_specs.pstring), "Mars")
        lim = 1.7;
        scatter(X2(1,1)/AU, X2(2,1)/AU, "red", 'filled');
    elseif  strcmp(string(P2_specs.pstring), "Venus")
        lim = 1.5;
        scatter(X2(1,1)/AU, X2(2,1)/AU, "green", 'filled');
    elseif strcmp(string(P2_specs.pstring), "Jupiter")
        lim = 5.5;
        scatter(X2(1,1)/AU, X2(2,1)/AU, MarkerFaceColor=[.5  0 .5], MarkerEdgeColor=[0 0 0]);
    elseif strcmp(string(P2_specs.pstring), "Mercury")
        lim = 1.5;
        scatter(X2(1,1)/AU, X2(2,1)/AU, MarkerFaceColor=[.7 .7 .7], MarkerEdgeColor=[0 0 0]);
    else
        lim = 1.7;
        scatter(X2(1,1)/AU, X2(2,1)/AU, "black", 'filled');
    end

    scatter(XSats_flat(1,:)/AU, XSats_flat(2,:)/AU, "black", "square");
    legend("Sun", "Earth", P2_specs.pstring, "Relay Satellites")
    hold off

    xlim([-lim, lim])
    ylim([-lim, lim])
    xlabel("x [AU]");
    ylabel("y [AU]");
    pbaspect([1 1 1])
    Nl = sum(bestGenome(2:length(bestGenome)/Nf:end));
    title("Current Best Constellation: " + int2str(Nl) + " Launches | " + int2str(sum(NSats)) + " Satellites Delivering an Average Bandwidth of " + num2str(-B) + " Mbps.")

end


% function [Out, XSats, NSats] = finish_with_fmincon(Xopt, X1, X2, etR, mu, Shuttle_specs, Comms_specs, Sat_specs, P1_specs, P2_specs, Nf, lb, ub, A, b)
% 
%     fun = @(X) wrapperFunc_RF_2P(X, X1, X2, etR, mu, Shuttle_specs, Comms_specs, Sat_specs, P1_specs, P2_specs, Nf);
% 
%     % Modify bounds to 
%     lb_intsRfixed 
% 
% 
% end