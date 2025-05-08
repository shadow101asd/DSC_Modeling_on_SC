function [A, b, lb, ub] = generateCircularConstraints(N, MIN_R, MAX_R, MAXSATS)
    A = zeros(3*N);
    b = zeros(3*N,1);
    
    % Sum of nums <= MAXSATS
    % A(1,1:N) = ones(N,1);
    % b(1,1) = MAXSATS;
    %Relative ordering of nums
    % if N >= 2
    %     for i = 1:(N-1)
    %         A(i+1,i) = -1;
    %         A(i+1,i+1) = 1;
    %     end
    % end
    
    lb = zeros(3*N,1);
    lb(N+1:2*N) = MIN_R*ones(N,1); % Each loop has to respect the min R
    lb(1) = 1; % We need at least 1 satellite in the simulation or things get weird

    ub = zeros(3*N,1);
    ub(1:N) = MAXSATS*ones(N,1); % Each loop has a max number of sats
    ub(N+1:2*N) = MAX_R*ones(N,1); % Each loop has to respect the max R
    ub(2*N+1:3*N) = 2*pi*ones(N,1);
end

