%% Perturbation Calcs

clear all

inputfile = "../Inputs/run008.mat";
datafile = "../Data/run008/D2/50.mat";
load(inputfile);
load(datafile);

muIndices = [1, 2, 3, 4, 5, 6];
% muIndices = 2;
mus = zeros(size(muIndices));

%% SPICE Data

ref = 'ECLIPJ2000';
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

% Mus
for i = 1:length(muIndices)
    mus(i) =  cspice_bodvcd(muIndices(i), 'GM', 10); % GM of the requested body with 10 significant digits
end

% Ephemerides
N = size(muIndices, 2);
nT = size(etR, 2);
XPs = zeros([6, nT, N]);
for i = 1:N
    idx = muIndices(i);
    XPs(:,:,i) = cspice_spkezr(int2str(idx), etR, ref, 'NONE', '10'); % Body ephemeris wrt to the central body in ref
end


%% Analysis

As = computeNetPerturbations(XSats_opt, XPs, mus, etR);
Ameans = mean(As, 2);
DVs = Ameans*(etR(end)- etR(1))*1e3 % m/s^2

%% Plotting

figure(1)
timeX = (etR-etR(1))/(24*3600);
plot(timeX, As)
set(gca, 'YScale', 'log')
legend(string(1:size(XSats_opt, 3)));