function y = Keplerian2Cartesian(a, e, i, Om, w, f0, mu)
    % Keplerian2Cartesian - Convert Keplerian elements to Cartesian state vectors
    %
    % Inputs (all 1×n):
    %   a   - semi-major axis
    %   e   - eccentricity or [3×n] eccentricity vectors
    %   i   - inclination [rad]
    %   Om  - RAAN [rad]
    %   w   - argument of periapsis [rad]
    %   f0  - true anomaly [rad]
    %   mu  - gravitational parameter (scalar or 1×n)
    % 
    % Alternative input: one Keplerian passed in
    %
    % Output:
    %   y   - 6×n matrix [r; v] in inertial frame

    if nargin == 2 % alternative case
        mu = e;
        [a,e,i,Om,w,f0] = unpackKeplerian(a);
    end

    % Handle eccentricity vector case
    if size(e, 1) == 3
        enorm = vecnorm(e);  % 1×n
    else
        enorm = e;
    end

    p = a.*(1-enorm.^2);
    n = size(a, 2);

    % Special Cases
    if norm(enorm) < 1e-6 * sqrt(n) % if circular
        w = 0*w;
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            %f0 = acos(r(1, :)./R); % lambda_true
            Om = Om * 0;
        else
            %f0 = acos(dot(N,r)./(Nnorm.*R)); % u
        end
    else
        if norm(i) < 1e-6 * sqrt(n) % if equatorial
            %f0 = acos(r(1, :)./R); % lambda_true
            Om = Om * 0;
        end
    end
    
    % Storing Variables for Computational Efficiency
    cosnu = cos(f0);
    sinnu = sin(f0);
    rootmup = sqrt(mu.*p.^-1);
    rPQW = [p.*cosnu./(1+ enorm.*cosnu); p.*sinnu./(1+ enorm.*cosnu); zeros(1,n)];
    vPQW = [-rootmup.*sinnu; rootmup.*(enorm + cosnu); zeros(1,n)];
    
    % Computing IJK/PQW
    IJKOverPQW1 = cos(Om).*cos(w) - sin(Om).*sin(w).*cos(i);
    IJKOverPQW2 = sin(Om).*cos(w) + cos(Om).*sin(w).*cos(i);
    IJKOverPQW3 = sin(w).*sin(i);
    IJKOverPQW4 = -cos(Om).*sin(w) - sin(Om).*cos(w).*cos(i);
    IJKOverPQW5 = -sin(Om).*sin(w) + cos(Om).*cos(w).*cos(i);
    IJKOverPQW6 = cos(w).*sin(i);
    IJKOverPQW7 = sin(Om).*sin(i);
    IJKOverPQW8 = -cos(Om).*sin(i);
    IJKOverPQW9 = cos(i);
    
    % for j = 1:n
    %     IJKOverPQW = [IJKOverPQW1(1, j), IJKOverPQW4(1, j), IJKOverPQW7(1, j);
    %                     IJKOverPQW2(1, j), IJKOverPQW5(1, j), IJKOverPQW8(1, j); 
    %                     IJKOverPQW3(1, j), IJKOverPQW6(1, j), IJKOverPQW9(1, j)];
    %     rPQW(:,j) = IJKOverPQW * rPQW(:,j);
    %     vPQW(:,j) = IJKOverPQW * vPQW(:,j);
    % end

    % Vectorized rotation of rPQW and vPQW into IJK
    rPQW_temp = rPQW;
    vPQW_temp = vPQW;

    rPQW(1,:) = IJKOverPQW1 .* rPQW_temp(1,:) + IJKOverPQW4 .* rPQW_temp(2,:) + IJKOverPQW7 .* rPQW_temp(3,:);
    rPQW(2,:) = IJKOverPQW2 .* rPQW_temp(1,:) + IJKOverPQW5 .* rPQW_temp(2,:) + IJKOverPQW8 .* rPQW_temp(3,:);
    rPQW(3,:) = IJKOverPQW3 .* rPQW_temp(1,:) + IJKOverPQW6 .* rPQW_temp(2,:) + IJKOverPQW9 .* rPQW_temp(3,:);

    vPQW(1,:) = IJKOverPQW1 .* vPQW_temp(1,:) + IJKOverPQW4 .* vPQW_temp(2,:) + IJKOverPQW7 .* vPQW_temp(3,:);
    vPQW(2,:) = IJKOverPQW2 .* vPQW_temp(1,:) + IJKOverPQW5 .* vPQW_temp(2,:) + IJKOverPQW8 .* vPQW_temp(3,:);
    vPQW(3,:) = IJKOverPQW3 .* vPQW_temp(1,:) + IJKOverPQW6 .* vPQW_temp(2,:) + IJKOverPQW9 .* vPQW_temp(3,:);


    y = [rPQW;vPQW];
end

