%% DSC Modular Setup File: Generate Ephemerides, and # of Satellites files for use in SC

clear all

%% SPICE-related files

path_to_generic_kernels = '/Users/jpenot/Documents/MATLAB/SPICE/naif.jpl.nasa.gov/pub/naif/generic_kernels';
path_to_mice            = '/Users/jpenot/Documents/MATLAB/SPICE/mice';
addpath(strcat(path_to_mice,'/src/mice'))
addpath(strcat(path_to_mice,'/lib'))

% Load the datafiles (kernels)
% (1) leap-seconds
cspice_furnsh([path_to_generic_kernels,'/lsk/naif0012.tls.pc']);

interval = 3600*24*3; % 3 days
ref = 'ECLIPJ2000';

% Dates
date0 = '2035 Jan 1 12:00:00 UTC'; % Simulation start
et0 = cspice_str2et(date0);
day = 24*3600; % one day in seconds
EM_synodic_T = 779.9 * day; % s
duration = EM_synodic_T * 3; % s
et1 = et0 + duration;
etR = et0:interval:et1; % Row vector of times between start and end date

[muSu, Xs] = generateSpiceEphemerides(ref, 10, etR, [3, 4]); % Earth and Mars

XEa = Xs(:,:,1);
XMa = Xs(:,:,2);


%% Saving

filepath = "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Inputs/";
filename = filepath + "run018";
save(filename, "muSu", "XEa", "XMa", "etR");

