function [new_arx, new_arfitness, new_arindex, new_lambda] = apply_extinction(arx, arfitness, arindex, p, min_lambda, k)
% Function keeps k best elements from arx array if k is set and randomly
% remove the rest
% If k is not set, fitness value is ignored and every element can be removed

if nargin < 6
  k = 0;
end

new_arx = zeros(size(arx));
new_arfitness = zeros(size(arfitness));
new_lambda = 0;

lambda = length(arfitness);
max_to_remove = lambda - min_lambda;
for i = 1:lambda
  if i < k || i > new_lambda + max_to_remove || rand < 1 - p 
    new_lambda = new_lambda + 1;
    new_arx(:, new_lambda) = arx(:,arindex(i));
    new_arfitness(new_lambda) = arfitness(i);
  end
end

new_arx = new_arx(:, 1:new_lambda);
new_arfitness = new_arfitness(1:new_lambda);
new_arindex = 1:new_lambda;

end
