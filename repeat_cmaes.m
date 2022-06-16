function [xmin, fitnessmin, out]=repeat_cmaes(fitness_function, dimensions, repetitions, use_best, extinction_type, lambda)
  % fitness_function - objective/fitness function 
  % dimensions - number of objective variables/problem dimension
  % extinction type - (0 - none, 1 - directed, 2 - random)
  if nargin < 6
    lambda = 310;
  end
  if nargin < 5
    extinction_type = 0;
  end

  rng(0) % Change 0 to another non-negative integer for different set of results or to 'shuffle' for random results each time
  seeds = randi(2^31, repetitions, 1);
  
  xmin = zeros(dimensions, 1);
  fitnessmin = 2^31;
  out.datx = [];
  out.arx = [];
  for i = 1:repetitions
    [temp_xmin, temp_fitnessmin, temp_out] = cmaes(fitness_function, dimensions, extinction_type, seeds(i), lambda);
    if use_best && temp_fitnessmin < fitnessmin
      xmin = temp_xmin;
      fitnessmin = temp_fitnessmin;
      out = temp_out;
    else
      % Update mean value in xmin
      xmin = ((i - 1) * xmin + temp_xmin) / i;
      % Update mean value in out.datx
      if height(out.datx) == 0
        out.datx = temp_out.datx;
      elseif height(temp_out.datx) == height(out.datx)
        out.datx = ((i-1) * out.datx + temp_out.datx) / i;
      elseif height(temp_out.datx) > height(out.datx)
        last = out.datx(end,:);
        diff = height(temp_out.datx) - height(out.datx);
        out.datx = [out.datx; repmat(last, diff , 1)];
        out.datx = ((i-1) * out.datx + temp_out.datx) / i;
      else
        last = temp_out.datx(end,:);
        diff = height(out.datx) - height(temp_out.datx);
        temp_out.datx = [temp_out.datx; repmat(last, diff , 1)];
        out.datx = ((i-1) * out.datx + temp_out.datx) / i;
      end
      % Update mean value in out.arx (same as above, but transposed and with different variable)
      if width(out.arx) == 0
        out.arx = temp_out.arx;
      elseif width(temp_out.arx) == width(out.arx)
        out.arx = ((i-1) * out.arx + temp_out.arx) / i;
      elseif width(temp_out.arx) > width(out.arx)
        last = out.arx(:, end);
        diff = width(temp_out.arx) - width(out.arx);
        out.arx = [out.arx repmat(last, 1, diff)];
        out.arx = ((i-1) * out.arx + temp_out.arx) / i;
      else
        last = temp_out.arx(:, end);
        diff = width(out.arx) - width(temp_out.arx);
        temp_out.arx = [temp_out.arx repmat(last, 1, diff)];
        out.arx = ((i-1) * out.arx + temp_out.arx) / i;
      end
    end
  end
  fitnessmin = fitness_function(xmin);
end