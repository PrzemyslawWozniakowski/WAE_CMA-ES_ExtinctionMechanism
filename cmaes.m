function [xmin, fitnessmin, out] = cmaes(fitness_function, dimensions, extinction_type, seed, lambda)
% fitness_function - objective/fitness function 
% dimensions - number of objective variables/problem dimension
% extinction type - (0 - none, 1 - directed, 2 - random)
  
  % --------------------  Initialization --------------------------------
  % Random numbers generator
  if nargin < 4
    seed = 'shuffle'; 
  end
  rng(seed)
  % User defined input parameters (need to be edited)
  xmean = rand(dimensions,1);    % objective variables initial point
  sigma = 0.5;          % coordinate wise standard deviation (step size)
  stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
  stopeval = 1e3*dimensions^2;   % stop after stopeval number of function evaluations

  % Strategy parameter setting: Selection
  if nargin < 5
    lambda = 310;
  end
  
  % Additional parameters for selection and adaptation
  [mu, weights, mueff, cc, cs, c1, cmu, damps] = update_params(lambda, dimensions);
  
  % Initialize dynamic (internal) strategy parameters and constants
  pc = zeros(dimensions,1); ps = zeros(dimensions,1);   % evolution paths for C and sigma
  B = eye(dimensions,dimensions);                       % B defines the coordinate system
  D = ones(dimensions,1);                      % diagonal D defines the scaling
  C = B * diag(D.^2) * B';            % covariance matrix C
  invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
  eigeneval = 0;                      % track update of B and D
  chiN=dimensions^0.5*(1-1/(4*dimensions)+1/(21*dimensions^2));  % expectation of
                                      %   ||N(0,I)|| == norm(randn(N,1))
  out.datx = [];  out.arx = []; % for plotting output

  % -------------------- Extinction settings --------------------------------
  c_extinction = 0.1;                   % threshold for difference between subsequent generations to call them stagnant
  p_extinction = 0.5;
  k_extinction = 0.2;
  count_stagnant = 0;                    % counter for currently stagnant generations
  extinction_trigger = 20;              % limit of stagnant generations which triggers extinction 
  min_lambda_fraction = 0.3;
  min_lambda = min_lambda_fraction * lambda;
  % -------------------- Generation Loop --------------------------------
  counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
  while counteval < stopeval

    % Generate and evaluate lambda offspring
    for k=1:lambda
      arx(:,k) = xmean + sigma * B * (D .* randn(dimensions,1)); % m + sig * Normal(0,C)
      arfitness(k) = feval(fitness_function, arx(:,k)); % objective function call
      counteval = counteval+1;
    end

    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness);  % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value

    % ----------------------- Extinction ----------------------------------
    if extinction_type ~= 0  && lambda > ceil(min_lambda) % If extinction should be considered
      if abs(norm(xmean-xold)) < c_extinction % If popuplation is stagnant, increase counter
        count_stagnant = count_stagnant + 1;
        if count_stagnant >= extinction_trigger
          old_lambda = lambda;
          if extinction_type == 1 % Targeted extinction
            [arx, arfitness, arindex, lambda] = apply_extinction(arx, arfitness, arindex, p_extinction, min_lambda, floor(k_extinction*lambda));
          elseif extinction_type == 2 % Random extinction
            [arx, arfitness, arindex, lambda] = apply_extinction(arx, arfitness, arindex, p_extinction, min_lambda);
          end
          if lambda ~= old_lambda % Update params for new population size after extinction
            [mu, weights, mueff, cc, cs, c1, cmu, damps] = update_params(lambda, dimensions); 
          end
        end
      else
        count_stagnant = 0; % Reset the counter if population is not stagnant
      end
    end
    % ---------------------------------------------------------------------

    % Cumulation: Update evolution paths
    ps = (1-cs) * ps ...
          + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/dimensions < 2 + 4/(dimensions+1);
    pc = (1-cc) * pc ...
          + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;

    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
    C = (1-c1-cmu) * C ...                   % regard old matrix
         + c1 * (pc * pc' ...                % plus rank one update
                 + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
         + cmu * artmp * diag(weights) * artmp'; % plus rank mu update

    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));

    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/dimensions/10  % to achieve O(N^2)
      eigeneval = counteval;
      C = triu(C) + triu(C,1)'; % enforce symmetry
      [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
      D = sqrt(diag(D));        % D contains standard deviations now
      invsqrtC = B * diag(D.^-1) * B';
    end

    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable
    if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
      break;
    end

    % Output
    out.datx = [out.datx; xmean'];
  end % while, end generation loop

  % ------------- Final Message and Plotting Figures --------------------
%   disp(['Oszacowań funkcji: ' num2str(counteval) ', Wynik: ' num2str(arfitness(1)) ', Typ wymarcia: ' num2str(extinction_type) ', Seed: ' num2str(seed)]);
  xmin = arx(:, arindex(1)); % Return best point of last iteration.
                             % Notice that xmean is expected to be even
                             % better.
  fitnessmin = arfitness(1); % Return fitness for bes point of last iteration.
  out.arx = arx;             % Return points from last generation.
end