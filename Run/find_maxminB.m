function [b, xsol] = find_maxminB(Bt, x0, numsats)

    % In this function, we assume each satellite has two antennas.

    NPs = size(Bt,1) - numsats - 1; % Number of planets (clients) other than Earth
    assert(NPs >= 2);
    
    % Set Up Optimization Problem

    

end