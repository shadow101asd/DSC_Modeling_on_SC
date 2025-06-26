function [Dtxs, Drxs, Ptxs, Nedges, eAs, Tsyss] = formSpecVectors2(NSats, Sat_specs, Earth_specs, Mars_specs, D_antennas)
    % Unpack relevant specs for link budget calculations, allowing for
    % different groups of satellites to have different antenna diameters
    
    % Incorporate different antenna diameters for different #s of
    % satellites
    DAs = reshape(repelem(D_antennas, NSats), [sum(NSats) 1]);
    Dtxs = [DAs; Earth_specs.Dtx; Mars_specs.Dtx];
    Drxs = [DAs; Earth_specs.Drx; Mars_specs.Drx];

    % Same as the original function
    Ptxs = [repelem(Sat_specs.Ptx, sum(NSats))'; Earth_specs.Ptx; Mars_specs.Ptx];
    Nedges = [repelem(Sat_specs.Nedges, sum(NSats))'; Earth_specs.Nedges; Mars_specs.Nedges];
    eAs = [repelem(Sat_specs.eA, sum(NSats))'; Earth_specs.eA; Mars_specs.eA];
    Tsyss = [repelem(Sat_specs.Tsys, sum(NSats))'; Earth_specs.Tsys; Mars_specs.Tsys];

end