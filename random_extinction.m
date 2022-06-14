function [new_arx, new_arfitness, new_arindex, new_lambda] = random_extinction(arx, arfitness, arindex, p, seed)
%Function randomly removes elements with probability p

if nargin < 5
    seed = 0;
end
rng(seed);

new_arx = zeros(1,length(arx));
new_arindex = zeros(1,length(arx));
new_arfitness = zeros(1,length(arx));
new_lambda = 0;

for i = 1:length(arx)
    if rand < 1 - p
        new_lambda = new_lambda + 1;
        new_arx(new_lambda) = arx(i);
        new_arindex = arindex(i);
        new_arfitness(new_lambda) = arfitness(arindex(i));
    end
end

new_arx = new_arx(1:new_lambda);
new_arindex = new_arindex(1:new_lambda);
new_arfitness = new_arfitness(1:new_lambda);

end