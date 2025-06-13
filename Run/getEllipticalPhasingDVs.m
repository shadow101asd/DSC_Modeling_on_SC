function [DVs, start_locs, arrival_locs] = getEllipticalPhasingDVs(a, e, mu, Dfs)

    % We need to consider phasing maneuvers from the periapsis and apoapsis
    % of our original elliptical orbit. Locs and dirs return arguments are
    % just for troubleshooting

    N = length(Dfs);

    rp = a*(1-e);
    ra = a*(1+e);
    
    % Use vis-viva
    vp = sqrt(mu*(2/rp - 1/a)); % orbital velocity at periapsis (in the original orbit)
    va = sqrt(mu*(2/ra - 1/a)); % orbital velocity at apoapsis  (in the original orbit)

    % Get phasing orbit options (to vectorize if necessary)

    a_tps = a*(1-Dfs).^(2/3);
    a_tas = a*(1+Dfs).^(2/3);

    et_p2p = 1 - a./a_tps*(1-e);
    et_p2a = 1 + a./a_tas*(1-e);
    et_a2p = 1 - a./a_tps*(1+e);
    et_a2a = 1 + a./a_tas*(1+e);
    
    % Use vis-viva again
    v_p2ps = sqrt(mu*(2/rp - 1./a_tps));
    v_p2as = sqrt(mu*(2/rp - 1./a_tas));
    v_a2ps = sqrt(mu*(2/ra - 1./a_tps));
    v_a2as = sqrt(mu*(2/ra - 1./a_tas));

    % Compute DVs

    DV_options = zeros([4, N]);
    DV_options(1,:) = 2*abs(vp - v_p2ps);
    DV_options(2,:) = 2*abs(vp - v_p2as);
    DV_options(3,:) = 2*abs(va - v_a2ps);
    DV_options(4,:) = 2*abs(va - v_a2as);

    % Find min DVs and report info

    [DVs, idxs] = min(DV_options);

    start_locs(N) = "periapsis";
    start_locs(idxs <= 2) = "periapsis";
    start_locs(idxs > 2) = "apoapsis";

    arrival_locs(N) = "periapsis";
    arrival_locs(idxs == 1) = "periapsis";
    arrival_locs(idxs == 3) = "periapsis";
    arrival_locs(idxs == 2) = "apoapsis";
    arrival_locs(idxs == 4) = "apoapsis";
end