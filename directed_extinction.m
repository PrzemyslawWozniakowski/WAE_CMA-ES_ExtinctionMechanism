function [new_arx, new_arfitness, new_arindex, new_lambda] = directed_extinction(arx, arfitness, arindex, k, p, seed)
%Function keeps k best elements from arx array

if nargin < 6
    seed = 0;
end

rng(seed)
k_best = arfitness(k);

new_arx = zeros(size(arx));
new_arindex = zeros(size(arindex));
new_arfitness = zeros(size(arfitness));
new_lambda = 0;

for i = 1:width(arx)
    if k_best > arfitness(arindex(i))
        if rand < 1 - p
            new_lambda = new_lambda + 1;
            new_arx(:, new_lambda) = arx(:,i);
            new_arindex = arindex(i);
            new_arfitness(new_lambda) = arfitness(arindex(i));
        end
    end
end

new_arx = new_arx(:, 1:new_lambda);
new_arindex = new_arindex(1:new_lambda);
new_arfitness = new_arfitness(1:new_lambda);

end

