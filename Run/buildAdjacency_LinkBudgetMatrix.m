function R = buildAdjacency_LinkBudgetMatrix(NSats, Comms_specs, Sat_specs, A_specs, B_specs)
    % To eventually generalize to more customers. A_lb is the datarate *
    % km^2 (just missing the FSPL essentially)

    c = 299792458; % speed of light, m/s
    lam = c/Comms_specs.f; % wavelength, m

    [Dtxs, Drxs, Ptxs, ~, eAs, Tsyss] = formSpecVectors(NSats, Sat_specs, A_specs, B_specs);

    % Compute Delivered Power at 1km distance
    Gtxs = (pi*Dtxs/lam).^2 .* eAs;
    Grxs = (pi*Drxs/lam).^2 .* eAs;

    FSPL_1km = (4*pi*1000/lam)^2; % FSPL over 1km. We then multiply our final matrix by 1/d^2 (d in km) to get the resulting bandwidth.

    Pr = Ptxs.* Gtxs * Grxs' / Comms_specs.Ls / Comms_specs.Lm / FSPL_1km; % W, Pr at 1 km, [Nsats + 2, Nsats + 2]

    % Noise power spectral density
    k_bolt = 1.380649e-23;
    N0s = k_bolt*Tsyss;

    B = Pr ./ N0s' / Comms_specs.SNR_desired; % N0s are relevant to the receiver side of the equation
    n = Comms_specs.encodingBpHz * Comms_specs.FEC; % bits/Hz encoding + FEC total spectral efficiency
    R = B*n;
end