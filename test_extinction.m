dimensions = 20;
% f = @(x)cec22_test_func(x, 6);
f = @(x)simple_functions.frastrigin10(x);
repetitions = 15;
lambda = 310;

% BASIC CMA-ES
[x, fitnessmin, out] = repeat_cmaes(f, dimensions, repetitions, false, 0, lambda, 0, 0);

disp(['Podstawowy CMA-ES, fmin = ' num2str(fitnessmin) ', argmin = ' mat2str(x')]);
figure(1);
hold off; 
plot(out.datx);
title("Krzywe zbieżności, funkcja CEC, podstawowy CMA-ES"); 
grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');

% figure(2);
% hold on;
% for n = 1:dimensions ecdf(out.arx(n, :)); end
% hold off;
% title("Dystrybuanta empiryczna, funkcja CEC, podstawowy CMA-ES"); 
% grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');

res_basic = fitnessmin;
res_random = zeros(5,6);
res_directed = zeros(5,6);

time_basic = out.elapsed;
time_random = zeros(5,6);
time_directed = zeros(5,6);

p_extinctions = [0.1 0.3 0.5 0.7 0.9];
extinction_triggers = 10:10:110;
k = 1;
for i = 1:length(p_extinctions)
  p_extinction = p_extinctions(i);
  for j = 1:length(extinction_triggers)
    extinction_trigger = extinction_triggers(j);
    % Random extinction
    [x, fitnessmin, out] = repeat_cmaes(f, dimensions, repetitions, false, 2, lambda, extinction_trigger, p_extinction);
%     disp(['Losowe wymieranie, K = ' num2str(extinction_trigger) ', p_e = ' num2str(p_extinction) ', fmin = ' num2str(fitnessmin) ', argmin = ' mat2str(x')]);
%     figure(3);
%     subplot(length(p_extinctions), length(extinction_triggers), k)
%     hold off; 
%     plot(out.datx);
%     title("Krzywe zbieżności, funkcja CEC, losowe wymieranie, K = " + extinction_trigger + "p_e = " + p_extinction); 
%     grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');
    res_random(i, j) = fitnessmin;
    time_random(i, j) = out.elapsed;
    
%     figure(4);
%     subplot(length(p_extinctions), length(extinction_triggers), k)
%     hold on;
%     for n = 1:dimensions 
%       ecdf(out.arx(n, :)); 
%     end
%     hold off;
%     title("Dystrybuanta empiryczna, funkcja CEC, losowe wymieranie, K = " + extinction_trigger + "p_e = " + p_extinction); 
%     grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');
%    
    % Targeted extinction
    [x, fitnessmin, out] = repeat_cmaes(f, dimensions, repetitions, false, 1, lambda, extinction_trigger, p_extinction);
%     disp(['Ukierunkowane wymieranie, K = ' num2str(extinction_trigger) ', p_e = ' num2str(p_extinction) ', fmin = ' num2str(fitnessmin) ', argmin = ' mat2str(x')]);
%     figure(5);
%     subplot(length(p_extinctions), length(extinction_triggers), k)
%     hold off; 
%     plot(out.datx);
%     title("Krzywe zbieżności, funkcja CEC, ukierunkowane wymieranie, K = " + extinction_trigger + "p_e = " + p_extinction); 
%     grid on; xlabel('Iteracje'); ylabel('Wartość x_i = argmin_i(f)');
    res_directed(i, j) = fitnessmin;
    time_directed(i, j) = out.elapsed;
%     
%     figure(6);
%     subplot(length(p_extinctions), length(extinction_triggers), k, 'align')
%     hold on;
%     for n = 1:dimensions 
%       ecdf(out.arx(n, :)); 
%     end
%     hold off;
%     title("Dystrybuanta empiryczna, funkcja Rastrigina, ukierunkowane wymieranie, K = " + extinction_trigger + "p_e = " + p_extinction); 
%     grid on; ylim([0 1.1]); xlabel('x_i'); ylabel('Wartość');


    k = k + 1;
  end
  disp(num2str(res_basic));
  disp(num2str(res_random));
  disp(num2str(res_directed));
% 
%   disp(time_basic);
%   disp(time_random);
%   disp(time_directed);
end
