function Kns = KnsfromK_2DFlower_Es(Ki,numsats,mu,W,M,E,emin,emax,distribution_name)
%KNSFROMK
    if nargin <= 5
        Kns = KnsfromK_2DFlower(Ki,numsats,mu,W,M);
        return
    end

    [a,e,i,Om,w,f0i] = unpackKeplerian(Ki);
    Kns(:,numsats) = zeros(6,1);

    if nargin == 6
        emin = e;
        emax = 1-e; % default behavior
    end

    if nargin <= 8
        distribution_name = "triangle";
    end

    % Assume orbital symmetry around the Sun for now
    for n = 1:numsats
        Kn = Ki;
        spacing_W = 2*pi/numsats * W;
        if e > 1e-4 % Orbit non-circular, change argument of periapsis w and eccentricity if applicable
            Kn(5) = w + spacing_W*(n-1);
            
            if E > 0
                switch distribution_name
                    case "triangle", distro_func = @distribute_eccentricities_triangle;
                    case "trapezoid",  distro_func = @distribute_eccentricities_trapezoid;
                    case "sine",  distro_func = @distribute_eccentricities_sine;
                    otherwise, error("Unknown distribution: %s", distribution_name);
                end
                Kn(2) = distro_func(E, numsats, n, emin, emax);
            end
        end

        % Regardless, stagger mean anomaly linearly out in time using
        % updateTrueAnomaly
        T = 2*pi*sqrt(a^3/mu); % Orbital Period
        dt = T + M*T*(n-1)/numsats;
        Kn(6) = updateTrueAnomaly(a,Kn(2),i,Om,w,f0i,mu,dt);

        Kns(:,n) = Kn;
    end

    function e = distribute_eccentricities_triangle(E, numsats, n, emin, emax)
    % DISTRIBUTE_ECCENTRICITIES
    % E       - number of roundtrips (min->max->min counts as 1)
    % numsats - total number of points
    % n       - (optional) 1-based index/indices to query; if empty/omitted, returns all
    % emin, emax - bounds for eccentricity
    
        idx = n(:).' - 1;            % convert to 0-based (accepts scalar or vector)
    
        R = emax - emin;
        denom = max(numsats-1, 1);       % avoid division by zero if numsats==1
    
        % Triangle wave in [0,1] using abs+mod, with E roundtrips over numsats samples
        x   = (2*E) * (idx / denom);     % phase that spans 0..2E across the sequence
        tri = 1 - abs(mod(x, 2) - 1);    % triangle wave in [0,1]
    
        % Map to [emin, emax]
        e = emin + R * tri;
    
        % Keep shape consistent with n input
        if nargin >= 3 && ~isempty(n)
            e = reshape(e, size(n));
        end
    end

    function e = distribute_eccentricities_sine(E, numsats, n, emin, emax)
    % DISTRIBUTE_ECCENTRICITIES_SINE
    % E       - number of roundtrips (min->max->min counts as 1)
    % numsats - total number of points
    % n       - (optional) 1-based index/indices to query; if empty/omitted, returns all
    % emin, emax - bounds for eccentricity
    
        if nargin < 3 || isempty(n)
            idx = 0:(numsats-1);          % all samples
        else
            idx = n(:).' - 1;             % convert to 0-based indices
        end
    
        R = emax - emin;
        denom = max(numsats-1, 1);
    
        % Phase goes from 0 .. 2*pi*E across numsats samples
        phase = 2*pi*E * (idx/denom);
    
        % Sine oscillates in [-1,1], scale to [0,1]
        s = (1 - cos(phase))/2;
    
        % Map to [emin, emax]
        e = emin + R * s;
    
        % Reshape if queried indices were passed
        if nargin >= 3 && ~isempty(n)
            e = reshape(e, size(n));
        end
    end
    
    function e = distribute_eccentricities_trapezoid(E, numsats, n, emin, emax, duty, edge)
    % DISTRIBUTE_ECCENTRICITIES_TRAPEZOID
    % Trapezoid-shaped eccentricity profile over E cycles across numsats samples.
    % Starts near emin at n=1, rises quickly to emax, holds, then falls quickly.
    %
    % Inputs:
    %   E       - number of roundtrips (min->max->min counts as 1)
    %   numsats - total number of samples
    %   n       - (optional) 1-based index/indices to query; if empty/omitted, returns all
    %   emin, emax - bounds for eccentricity
    %   duty    - (optional) fraction of each cycle held at emax (default 0.85)
    %   edge    - (optional) fraction per ramp (rise and fall each) (default 0.075)
    %
    % Constraints: 0 <= edge <= 0.5, 0 <= duty <= 1 - 2*edge
    % The remaining fraction (low_hold) = 1 - (2*edge + duty) is the time at emin.
    
        if nargin < 6 || isempty(duty), duty = 0.85; end
        if nargin < 7 || isempty(edge), edge = 0.075; end
    
        if edge < 0 || edge > 0.5
            error('edge must be in [0, 0.5].');
        end
        if duty < 0 || duty > 1 - 2*edge
            error('duty must satisfy 0 <= duty <= 1 - 2*edge.');
        end
    
        % Indices: all or a requested subset
        if nargin < 3 || isempty(n)
            idx = 0:(numsats-1);            % all samples
        else
            idx = n(:).' - 1;               % 0-based indices
        end
    
        R = emax - emin;
        denom = max(numsats-1, 1);
    
        % Phase in [0,1) repeating E times across the sequence
        phi = mod(E * (idx / denom), 1);
    
        % Trapezoid parameters
        low_hold = 1 - (2*edge + duty);     % time at emin (can be zero)
        t0 = low_hold;
        t1 = t0 + edge;                     % end of rising edge
        t2 = t1 + duty;                     % end of high plateau
        t3 = t2 + edge;                     % end of falling edge (=1 when low_hold==0)
    
        % Build trapezoid in [0,1]
        s = zeros(size(phi));               % start at low
        rise = (phi >= t0) & (phi < t1);
        high = (phi >= t1) & (phi < t2);
        fall = (phi >= t2) & (phi < t3);
    
        s(rise) = (phi(rise) - t0) / max(edge, eps);
        s(high) = 1;
        s(fall) = 1 - (phi(fall) - t2) / max(edge, eps);
        % (Remaining region, if any, is low plateau at 0)
    
        % Map to [emin, emax]
        e = emin + R * s;
    
        % Preserve shape if user queried specific indices
        if nargin >= 3 && ~isempty(n)
            e = reshape(e, size(n));
        end
    end
end