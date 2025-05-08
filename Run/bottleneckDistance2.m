function Dmax = bottleneckDistance2(graph,edgepath)
%BOTTLENECKDISTANCE Summary of this function goes here
%   Obtain the bottleneck distance (maximum weight edge) along a path in a
%   graph.

    % Extract the weights and take the maximum
    Dmax = max(graph.Edges.Weight(edgepath));
end

