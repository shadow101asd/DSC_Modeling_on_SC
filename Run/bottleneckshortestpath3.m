function [path, D, edgepath, D_Bottleneck] = bottleneckshortestpath3(g, nodename1, nodename2, reconstruct_path_bool)
%BOTTLENECKSHORTESTPATH
    T = minspantree(g, Method="sparse", Root=nodename1);
    [path1, D1, edgepath_t] = mst_path(T,nodename1,nodename2);


    D_Bottleneck_temp = bottleneckDistance2(T,edgepath_t);

    if reconstruct_path_bool
        % Return and check for more direct, roughly equivalent bottleneck paths
        margin = 1.02;
        DMax = margin*D_Bottleneck_temp;
        
        edges2remove = find(g.Edges.Weight > DMax);
        graph2 = rmedge(g, edges2remove);
    
        [path,D,edgepath] = shortestpath(graph2,nodename1,nodename2);
        D_Bottleneck = bottleneckDistance2(graph2, edgepath); % Edgepath is wrong here!
    else
        D_Bottleneck = D_Bottleneck_temp;
        path = string(g.Nodes.Name(path1)); % still report non-optimized path
        D = D1;
        edgepath = [];
    end
end

