dimensions = 13;

for extinction_type = 0:2
  [x, out] = purecmaes(@simple_functions.fsphere, dimensions, extinction_type);

  % Convergence curve
  figure(1);
  subplot(2,2,extinction_type + 1, 'align')
  hold off; 
  plot(out.datx);
  title("Krzywe zbieżności, typ wymarcia = " + extinction_type); 
  grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');
  
  % Empirical cumulative distribution function plot
  figure(2);
  subplot(2,2,extinction_type + 1)
  hold on;
  for n = 1:dimensions 
    ecdf(out.arx(n, :)); 
  end
  hold off;
  title("Dystrybuanta empiryczna dla ostatniej iteracji, typ wymarcia = " + extinction_type); 
  grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');
end

% purecmaes(@simple_functions.fsphere, 1);
% purecmaes(@simple_functions.fsphere, 2);

% Example for calling inline function
% purecmaes(@(x)(sum(x.^2)))

% Example for calling some function from variable
% test_fun = @(x) (sum(x.^2));
% purecmaes(test_fun)
