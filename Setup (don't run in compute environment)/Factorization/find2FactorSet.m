function array = find2FactorSet(int)
%FIND2FACTORSET Summary of this function goes here
%   Detailed explanation goes here
array = [];
fs = [1, factor(int)];

N = length(fs);

for k = 1:N
    C = unique(nchoosek(fs, k), "rows", "stable");
    
    for i = 1:size(C, 1)
        p1 = prod(C(i,:));
        p2 = int/p1;
        array = cat(1, array, [p1,p2]);
    end
end

array = unique(array, 'sorted', 'rows');
end
