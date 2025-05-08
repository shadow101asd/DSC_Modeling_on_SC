%% Generate 3-factor sets
clear all
MIN_NSATS = 1;
MAX_NSATS = 9600;

% Load memotable
memotable = "../../Inputs/MemoTables/3Factors";
load(memotable, "ThreeFs");

tic
for i = MIN_NSATS:MAX_NSATS
    if ~isfield(ThreeFs, "trip"+int2str(i))
        % Then there is no corresponding entry in the memo table, and we
        % should generate one
        t = find3FactorSet(i);
        ThreeFs.("trip"+int2str(i)) = t;
    end
end

save(memotable, "ThreeFs");

toc