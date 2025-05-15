function flag = bfs_reachable_At(At, s, t)
% Fastest possible BFS using adjacency matrix At with pointer queue
% At must be filtered to include only valid edges (nonzero entries)

    N = size(At, 1);
    visited = false(N, 1);
    queue = zeros(N, 1);    % preallocate queue to max possible size
    head = 1;
    tail = 1;

    queue(tail) = s;
    visited(s) = true;

    while head <= tail
        u = queue(head); head = head + 1;
        if u == t
            flag = true;
            return
        end
        % Get neighbors directly from sparse adjacency matrix row
        nbrs = find(At(u, :) ~= 0);
        for v = nbrs
            if ~visited(v)
                tail = tail + 1;
                queue(tail) = v;
                visited(v) = true;
            end
        end
    end

    flag = false;
end