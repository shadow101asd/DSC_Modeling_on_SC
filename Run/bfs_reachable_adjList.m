function flag = bfs_reachable_adjList(adjList, s, t)
% BFS connectivity test on prebuilt adjList

    visited = false(length(adjList), 1);
    queue = s;
    visited(s) = true;

    while ~isempty(queue)
        u = queue(1); queue(1) = [];
        if u == t
            flag = true;
            return
        end
        for v = adjList{u}
            if ~visited(v)
                visited(v) = true;
                queue(end+1) = v;
            end
        end
    end

    flag = false;
end
