function p = elliptical_orbit_perimeter(a,e)
%ELLIPTICAL_ORBIT_PERIMETER
%Numerically computes the perimeter/length of a single elliptical orbit of
% semimajor axis a and eccentricity e
% a, p are intended to be in [km], although in practice p just takes the
% units of a
% e is unitless

b = a*sqrt(1-e^2);

% Derivative of the parametric equations of the ellipse

dxdt = @(t) -a * sin(t);
dydt = @(t) b * cos(t);

% Calculate perimeter using numerical integration

p = integral(@(t) sqrt(dxdt(t).^2 + dydt(t).^2), 0, pi, 'ArrayValued', true) * 2;

end