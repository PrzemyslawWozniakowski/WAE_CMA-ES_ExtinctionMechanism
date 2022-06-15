function [new_arx, new_arfitness, new_arindex, new_lambda] = extinction(arx, arfitness, arindex, p, seed, k)
% Function keeps k best elements from arx array if k is set
% Otherwise, fitness value is ignored and every element can be removed

if nargin < 6 || k == length(arx)
  ignore_fitness = true;
else
  ignore_fitness = false;
  k_best = arfitness(k);
end

if nargin < 5
  seed = 0;
end

rng(seed)

new_arx = zeros(size(arx));
new_arindex = zeros(size(arindex));
new_arfitness = zeros(size(arfitness));
new_lambda = 0;

for i = 1:width(arx)
    if ignore_fitness || k_best > arfitness(arindex(i))
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

