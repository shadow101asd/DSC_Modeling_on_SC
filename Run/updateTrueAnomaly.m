function f0new = updateTrueAnomaly(a, e, i, Om, w, f0, mu, dt)
%UPDATETRUEANOMALY Summary of this function goes here
%   Detailed explanation goes here

%#codegen

% Compute current eccentric anomaly

Ecur = 2*atan(sqrt((1-e)/(1+e))*tan(f0/2));

% Compute Current Mean Anomaly from Kepler's Formula

Mecur = Ecur - e*sin(Ecur);

% Update Mean Anomaly

n = sqrt(mu/a^3); % Mean motion
Meupdated = Mecur + n*dt;

% Compute Updated Eccentric Anomaly from Meupdated
% Numerical method:

% fun = @(E) E - e*sin(E) - Meupdated;
% Eupdated = fzero(fun, Ecur); % Ecur used as initial guess


% Use Newton-Raphson: Runs faster than fzero afaict

tol = 1e-5;
max_iter = 15;

% Initial guess (Mupdated is good when e is small)
Eupdated = Meupdated;

for k = 1:max_iter
    f = Eupdated - e*sin(Eupdated) - Meupdated;
    fp = 1 - e*cos(Eupdated);
    dE = -f / fp;
    Eupdated = Eupdated + dE;

    % Convergence criteria
    if abs(dE) < tol
        % Compute Updated True Anomaly from Eupdated
        f0new = atan2((sin(Eupdated)*sqrt(1-e^2)),(cos(Eupdated)-e));
        return;
    end
end

% Warn if not converged, but still report result
f0new = atan2((sin(Eupdated)*sqrt(1-e^2)),(cos(Eupdated)-e));
warning('Newton method did not converge. Residual: %g', abs(f));
end
