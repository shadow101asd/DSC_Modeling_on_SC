function array = find3FactorSet(int)
%FIND3FACTORSET Summary of this function goes here
%   Detailed explanation goes here
array = [];
fs = [1, 1, factor(int)];
N = length(fs);

for k = 1:N
    C = unique(nchoosek(fs, k), "rows", "stable");
    
    for i = 1:size(C, 1)
        p1 = prod(C(i,:));
        p23 = find2FactorSet(int/p1);
        for j = 1:size(p23,1)
            r = p23(j,:);
            p2 = r(1);
            p3 = r(2);
            array = cat(1, array, [p1,p2,p3]);
        end
    end
end

array = unique(array, 'sorted', 'rows');
end

