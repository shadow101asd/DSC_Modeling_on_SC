function [] = DSC4_circ_SC_run(idx, run_idx)

% Run details
% Circular loops optimization
% SATNUMS = 5:1:200

% Constants

AU = 1.496e8; % 1 AU in km

% Load ephemerides and other pregenerated data
load("../Inputs/run"+run_idx+".mat", "muSu", "XEa", "XMa", "SATNUMS", "etR");

MAXSATS = SATNUMS(idx)
% Num timesteps
[~,nT] = size(XEa);

% Genetic Algorithm

% Constraints!

MIN_R = 0.5; % AU
MAX_R = 1.5; % AU

MAX_LOOPS = 5; % Has to be a positive integer
intcon = 1:MAX_LOOPS;

[A, b, lb, ub] = generateCircularConstraints(MAX_LOOPS, MIN_R, MAX_R, MAXSATS);

options = optimoptions('ga','Display','iter','FunctionTolerance', 1e-3, ...
         'MaxStallGenerations', 7, 'CreationFcn', 'gacreationuniformint',...
         'PopulationSize', 1000);

% Running the GA:
numvars = 3*MAX_LOOPS;
    
[X_opt, fval, EXIT_FLAG, OUTPUT] = ga(@(X) wrapperFuncCirc(XEa,XMa,X,etR,muSu),numvars, A, b, [],[],lb,ub,@(X)nonlcon(X, MAXSATS),intcon,options);

X_opt

% Saving results

filename = "../Data/run"+run_idx+"/" + int2str(MAXSATS);
save(filename, "X_opt", "fval", "EXIT_FLAG", "OUTPUT")

% Nested Functions

function Out = wrapperFuncCirc(X1,X2,X,etR,muSu)
    
    [nums, as, f0s] = unpackVarsCirc(X);
    Kns = KnsCircs(nums, as, f0s, muSu);
    [~, nT] = size(etR);
    
    numsats = sum(nums);
    XSats(:,:,numsats) = zeros(6,nT);
    for n = 1:numsats
        Kn = Kns(:,n);
        XSats(:,:,n) = propagateFromKeplerians(Kn,muSu,etR);
    end

    Out = bestLinkBudget(X1,X2,XSats);
end

function Kns = KnsCircs(nums, as, f0s, mu)
    N = length(nums);
    i = 1;
    Ki = [as(i), 0, 0, 0, 0, f0s(i)];
    Kns = Circular_Kns(Ki, nums(i), mu);
    for i = 2:N
        if nums(i) > 0
            Ki = [as(i), 0, 0, 0, 0, f0s(i)];
            Kns = cat(2, Kns, Circular_Kns(Ki, nums(i), mu));
        end
    end
end

function [nums, as, f0s] = unpackVarsCirc(X)
    AU = 1.496e8; % 1 AU in km
    % Unpack variables and convert as to km
    N = length(X)/3;
    nums = X(1:N);
    as = X(N+1:2*N)*AU;
    f0s = X(2*N+1:3*N);
end

function metric = bestLinkBudget(X1,X2,XSats)
    AU = 1.496e8; % 1 AU in km 
    % TLDR: What's the metric that we're optimizing over?
    [~,~,~,~,Dequivs] = networkAnalysis(X1,X2,"Earth","Mars",XSats);
    metric = mean(Dequivs)/AU;
end

function [c,ceq] = nonlcon(X,MAXSATS)
% This is technically a linear constraint, but I'll keep it here for now
% for clarity.
    [nums, ~, ~] = unpackVarsCirc(X);
    numsats = sum(nums);
    ceq = []; % No equality constraints for MINLPs
    c = numsats - MAXSATS;
end

end

