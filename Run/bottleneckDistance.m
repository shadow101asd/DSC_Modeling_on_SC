function Dmax = bottleneckDistance(graph,path)
%BOTTLENECKDISTANCE Summary of this function goes here
%   Obtain the bottleneck distance (maximum weight edge) along a path in a
%   graph.

    % Get edge indices along the path
    edgeIdx = findedge(graph, path(1:end-1), path(2:end));

    % Extract the weights and take the maximum
    if all(edgeIdx, 'all')
        Dmax = max(graph.Edges.Weight(edgeIdx));
    else
        Dmax = inf; % Path is infeasible in the passed-in graph
    end
end

