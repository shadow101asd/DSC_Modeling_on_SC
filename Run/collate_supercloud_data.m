function [] = collate_supercloud_data(run_idx)
% run_idx must be of the form "XXX", including leading zeros

myDir = "../Data/run" + run_idx + "/";


myFiles = dir(fullfile(myDir,'*.mat')); % gets all mat files in struct


N = length(myFiles);
satnums = zeros([1 N]);
fvals = zeros([1 N]);
exitflags = zeros([1 N]);
outputs = [];
X_opts = [];

for k = 1:N
    % Load
    baseFileName = myFiles(k).name;
    fullFileName = fullfile(myDir, baseFileName);
    load(fullFileName, 'X_opt', 'fval', 'EXIT_FLAG', 'OUTPUT');

    % Append
    str = split(baseFileName, '.');
    satnums(k) = str2double(str{1});

    X_opts = cat(3, X_opts, X_opt);
    fvals(k) = fval;
    exitflags(k) = EXIT_FLAG;
    outputs = cat(1, outputs, OUTPUT);

    % Clear to prevent errors
    clear X_opt fval EXIT_FLAG OUTPUT
end

% Reorder in increasing satnum order

[satnums, I] = sort(satnums);
X_opts = X_opts(:,:, I);
fvals = fvals(I);
exitflags = exitflags(I);
outputs = outputs(I);


%% Save

save_location = "../PP_Data/";
filename = "run" + run_idx;
save(save_location+filename, 'satnums', 'X_opts', 'fvals', 'exitflags', 'outputs')

end