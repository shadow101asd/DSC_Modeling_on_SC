function Dmax = bottleneckDistance(graph,path)
%BOTTLENECKDISTANCE Summary of this function goes here
%   Obtain the bottleneck distance (maximum weight edge) along a path in a
%   graph.

    % Get edge indices along the path
    edgeIdx = findedge(graph, path(1:end-1), path(2:end));

    % Extract the weights and take the maximum
    Dmax = max(graph.Edges.Weight(edgeIdx));
end

