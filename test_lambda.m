dimensions = 15;

for lambda = 10:50:410
  [x, fitnessmin, out] = repeat_cmaes(@simple_functions.frosenbrock, dimensions, 15, false, 0, lambda);

  % Convergence curve
  disp(['Podstawowy CMA-ES, lambda = ' num2str(lambda)  ', fmin = ' num2str(fitnessmin) ', argmin = ' mat2str(x')]);
  figure(1);
  subplot(3,3,floor(lambda/50) + 1, 'align')
  hold off; 
  plot(out.datx);
  title("Krzywe zbieżności, funkcja Rosenbrocka, lambda = " + lambda); 
  grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');
  
  % Empirical cumulative distribution function plot
  figure(2);
  subplot(3,3,floor(lambda/50) + 1)
  hold on;
  for n = 1:dimensions 
    ecdf(out.arx(n, :)); 
  end
  hold off;
  title("Dystrybuanta empiryczna, funkcja Rosenbrocka, lambda = " + lambda); 
  grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');
end
