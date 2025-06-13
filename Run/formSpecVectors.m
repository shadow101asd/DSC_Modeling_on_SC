function [Dtxs, Drxs, Ptxs, Nedges, eAs, Tsyss] = formSpecVectors(NSats, Sat_specs, Earth_specs, Mars_specs)
    % Unpack relevant specs for link budget calculations

    Dtxs = [repmat(Sat_specs.Dtx, [NSats, 1]); Earth_specs.Dtx; Mars_specs.Dtx];
    Drxs = [repmat(Sat_specs.Drx, [NSats, 1]); Earth_specs.Drx; Mars_specs.Drx];
    Ptxs = [repmat(Sat_specs.Ptx, [NSats, 1]); Earth_specs.Ptx; Mars_specs.Ptx];
    Nedges = [repmat(Sat_specs.Nedges, [NSats, 1]); Earth_specs.Nedges; Mars_specs.Nedges];
    eAs = [repmat(Sat_specs.eA, [NSats, 1]); Earth_specs.eA; Mars_specs.eA];
    Tsyss = [repmat(Sat_specs.Tsys, [NSats, 1]); Earth_specs.Tsys; Mars_specs.Tsys];

end