function [xmin, fitnessmin, out]=repeat_cmaes(fitness_function, dimensions, repetitions, use_best, extinction_type, lambda, extinction_trigger, p_extinction)
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
  out.elapsed = 0;
  for i = 1:repetitions
    [temp_xmin, temp_fitnessmin, temp_out] = cmaes(fitness_function, dimensions, extinction_type, seeds(i), lambda, extinction_trigger, p_extinction);
    out.elapsed = out.elapsed + temp_out.elapsed;
    if use_best && temp_fitnessmin < fitnessmin
      xmin = temp_xmin;
      fitnessmin = temp_fitnessmin;
      out = temp_out;
    elseif ~use_best
      % Update mean value in xmin
      xmin = xmin + temp_xmin;
      % Update mean value in out.datx
      if height(out.datx) == 0
        out.datx = temp_out.datx;
      elseif height(temp_out.datx) == height(out.datx)
        out.datx = out.datx + temp_out.datx;
      elseif height(temp_out.datx) > height(out.datx)
        last = out.datx(end,:);
        diff = height(temp_out.datx) - height(out.datx);
        out.datx = [out.datx; repmat(last, diff , 1)];
        out.datx = out.datx + temp_out.datx;
      else
        last = temp_out.datx(end,:);
        diff = height(out.datx) - height(temp_out.datx);
        temp_out.datx = [temp_out.datx; repmat(last, diff , 1)];
        out.datx = out.datx + temp_out.datx;
      end
      % Update mean value in out.arx (same as above, but transposed and with different variable)
      if width(out.arx) == 0
        out.arx = temp_out.arx;
      elseif width(temp_out.arx) == width(out.arx)
        out.arx = out.arx + temp_out.arx;
      elseif width(temp_out.arx) > width(out.arx)
        last = out.arx(:, end);
        diff = width(temp_out.arx) - width(out.arx);
        out.arx = [out.arx repmat(last, 1, diff)];
        out.arx = out.arx + temp_out.arx;
      else
        last = temp_out.arx(:, end);
        diff = width(out.arx) - width(temp_out.arx);
        temp_out.arx = [temp_out.arx repmat(last, 1, diff)];
        out.arx = out.arx + temp_out.arx;
      end
    end
  end

  out.elapsed = out.elapsed / repetitions;
  if ~use_best
    xmin = xmin / repetitions;
    out.arx = out.arx / repetitions;
    out.datx = out.datx / repetitions;
    fitnessmin = fitness_function(xmin);
  end
end