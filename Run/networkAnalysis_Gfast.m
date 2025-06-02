function [Dequivs, D_bottleneck_avg] = networkAnalysis_Gfast(X1,X2,~,~,XSats,alpha)
    AU = 149600000; % AU in km
    
    % Set up paths: This needs to be fixed to run on SC or elsehwere than this Mac
    gurobi_file_path = horzcat('/Library/gurobi1201/macos_universal2', '/matlab/gurobi_setup.m');
    gurobi_add_path = horzcat('/Library/gurobi1201/macos_universal2', '/examples/matlab/');
    
    % Setup gurobi
    run(gurobi_file_path)
    addpath(gurobi_add_path); % temporarily overwrites MATLAB's 'intlinprog.m'

    % Create Adjacency Matrices
    Xs = collateXs(X1,X2,XSats)/AU;
    D12_i = distanceBetweenXs(X1(:,1),X2(:,1))/AU;
    A = createAdjacencyMatrix2(Xs,@distanceBetweenXs);
    nT = size(Xs,2);

    Dequivs = zeros(nT,1);
    
    t = 1;
    while t <= nT
        At = A(:,:,t);
        if t == 1
            At(At>D12_i) = 0.0;
        else
            At(At>alpha*Dequivs(t-1)) = 0.0;
        end
        
        % Compute bottleneck shortest path
        [D_bottleneck, exitflag] = bottleneckshortestpath_GUROBI(At);

        % Catch errors
        if exitflag ~= 1 
            % Adjust alpha and try again
            alpha = alpha + 0.1;
            continue
        else % We're good to go
            Dequivs(t) = D_bottleneck;
            t = t + 1; % Move on to next loop
        end
        
    end

    D_bottleneck_avg = mean(Dequivs);
end