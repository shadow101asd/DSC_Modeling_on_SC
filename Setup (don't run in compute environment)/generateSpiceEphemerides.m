function [mu, Xs] = generateSpiceEphemerides(ref, muID, etR, ephemeridesArrayofIndices)
%GENERATESPICEEPHEMERIDES Summary of this function goes here
%   Detailed explanation goes here

path_to_generic_kernels = '/Users/jpenot/Documents/MATLAB/SPICE/naif.jpl.nasa.gov/pub/naif/generic_kernels';
path_to_mice            = '/Users/jpenot/Documents/MATLAB/SPICE/mice';
addpath(strcat(path_to_mice,'/src/mice'))
addpath(strcat(path_to_mice,'/lib'))

% Load the datafiles (kernels)
% (1) leap-seconds
cspice_furnsh([path_to_generic_kernels,'/lsk/naif0012.tls.pc']);
% (2) planets
cspice_furnsh([path_to_generic_kernels,'/spk/planets/de430.bsp']);
% (3) gravity constants
cspice_furnsh([path_to_generic_kernels,'/pck/gm_de431.tpc']);
% (4) planetary constant - you can also open this file with a text editor
% to read its contentmex -setup
cspice_furnsh([path_to_generic_kernels,'/pck/pck00010.tpc']);

% Pull data through SPICE

% Mu
mu = cspice_bodvcd(muID, 'GM', 10); % GM of the requested body with 10 significant digits

% Ephemerides
[~, N] = size(ephemeridesArrayofIndices);
[~, nT] = size(etR);
Xs = zeros([6, nT, N]);
for i = 1:N
    idx = ephemeridesArrayofIndices(i);
    Xs(:,:,i) = cspice_spkezr(int2str(idx), etR, ref, 'NONE', '10'); % Body ephemeris wrt to the central body in ref
end
end

