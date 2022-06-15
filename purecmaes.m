function xmin=purecmaes
  % --------------------  Initialization --------------------------------
  % User defined input parameters (need to be edited)
  strfitnessfct = 'frosenbrock';  % name of objective/fitness function
  N = 13;               % number of objective variables/problem dimension
  xmean = rand(N,1);    % objective variables initial point
  sigma = 0.5;          % coordinate wise standard deviation (step size)
  stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
  stopeval = 1e3*N^2;   % stop after stopeval number of function evaluations

  % Strategy parameter setting: Selection
  lambda = 4+floor(3*log(N));  % population size, offspring number
  
  % Additional parameters for selection and adaptation
  [mu, weights, mueff, cc, cs, c1, cmu, damps] = update_params(lambda, N);
  
  % Initialize dynamic (internal) strategy parameters and constants
  pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
  B = eye(N,N);                       % B defines the coordinate system
  D = ones(N,1);                      % diagonal D defines the scaling
  C = B * diag(D.^2) * B';            % covariance matrix C
  invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
  eigeneval = 0;                      % track update of B and D
  chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of
                                      %   ||N(0,I)|| == norm(randn(N,1))
  out.datx = [];  % for plotting output

  % -------------------- Extinction settings --------------------------------
  c_extinction = 0.1;                   % threshold for difference between subsequent generations to call them stagnant
  extinction_type = 2;                  % extinction type (0 - none, 1 - directed, 2 - random)
  p_extinction = 0.9;
  k_extinction = 0.75;
  count_stagnant = 0;                    % counter for currently stagnant generations
  extinction_trigger = 100;              % limit of stagnant generations which triggers extinction 
  min_lambda_fraction = 0.3;
  min_lambda = min_lambda_fraction * lambda;
  % -------------------- Generation Loop --------------------------------
  counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
  while counteval < stopeval

    % Generate and evaluate lambda offspring
    for k=1:lambda
      arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
      arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call
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
            [arx, arfitness, arindex, lambda] = extinction(arx, arfitness, arindex, p_extinction, min_lambda, 0, floor(k_extinction*lambda));
          elseif extinction_type == 2 % Random extinction
            [arx, arfitness, arindex, lambda] = extinction(arx, arfitness, arindex, p_extinction, min_lambda, 0);
          end
          if lambda ~= old_lambda % Update params for new population size after extinction
            [mu, weights, mueff, cc, cs, c1, cmu, damps] = update_params(lambda, N); 
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
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
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
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
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
  disp(['Ostatnia iteracja: ' num2str(counteval) ', Wynik: ' num2str(arfitness(1))]);
  xmin = arx(:, arindex(1)); % Return best point of last iteration.
                             % Notice that xmean is expected to be even
                             % better.

  % Convergence curve
  figure(1); hold off; 
  plot(out.datx);
  title('Krzywe zbieżności'); 
  grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');

  % Empirical cumulative distribution function plot
  figure(2); hold on;
  for n = 1:N 
    ecdf(arx(n, :)); 
  end
  hold off;
  title('Dystrybuanta empiryczna dla ostatniej iteracji'); 
  grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');

% ---------------------------------------------------------------
% Functions for test purposes
function f=frosenbrock(x)
  if size(x,1) < 2 error('dimension must be greater one'); end
  f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);

function f=fsphere(x)
  f=sum(x.^2);

function f=fssphere(x)
  f=sqrt(sum(x.^2));

function f=fschwefel(x)
  f = 0;
  for i = 1:size(x,1)
    f = f+sum(x(1:i))^2;
  end

function f=fcigar(x)
  f = x(1)^2 + 1e6*sum(x(2:end).^2);

function f=fcigtab(x)
  f = x(1)^2 + 1e8*x(end)^2 + 1e4*sum(x(2:(end-1)).^2);

function f=ftablet(x)
  f = 1e6*x(1)^2 + sum(x(2:end).^2);

function f=felli(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  f=1e6.^((0:N-1)/(N-1)) * x.^2;

function f=felli100(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  f=1e4.^((0:N-1)/(N-1)) * x.^2;

function f=fplane(x)
  f=x(1);

function f=ftwoaxes(x)
  f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);

function f=fparabR(x)
  f = -x(1) + 100*sum(x(2:end).^2);

function f=fsharpR(x)
  f = -x(1) + 100*norm(x(2:end));

function f=fdiffpow(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));

function f=frastrigin10(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  scale=10.^((0:N-1)'/(N-1));
  f = 10*size(x,1) + sum((scale.*x).^2 - 10*cos(2*pi*(scale.*x)));

function f=frand(x)
  f=rand;