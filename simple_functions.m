classdef simple_functions
  methods (Static)
    % Functions for test purposes
    function f=frosenbrock(x)
      if size(x,1) < 2 
        error('dimension must be greater one'); 
      end
      f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
    end
    
    function f=fsphere(x)
      f=sum(x.^2);
    end
    
    function f=fssphere(x)
      f=sqrt(sum(x.^2));
    end
    
    function f=fschwefel(x)
      f = 0;
      for i = 1:size(x,1)
        f = f+sum(x(1:i))^2;
      end
    end
    
    function f=fcigar(x)
      f = x(1)^2 + 1e6*sum(x(2:end).^2);
    end
    
    function f=fcigtab(x)
      f = x(1)^2 + 1e8*x(end)^2 + 1e4*sum(x(2:(end-1)).^2);
    end
    
    function f=ftablet(x)
      f = 1e6*x(1)^2 + sum(x(2:end).^2);
    end
    
    function f=felli(x)
      N = size(x,1); 
      if N < 2 
        error('dimension must be greater one'); 
      end
      f=1e6.^((0:N-1)/(N-1)) * x.^2;
    end
    
    function f=felli100(x)
      N = size(x,1); if N < 2 error('dimension must be greater one'); end
      f=1e4.^((0:N-1)/(N-1)) * x.^2;
    end
    
    function f=fplane(x)
      f=x(1);
    end
    
    function f=ftwoaxes(x)
      f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);
    end
    
    function f=fparabR(x)
      f = -x(1) + 100*sum(x(2:end).^2);
    end
    
    function f=fsharpR(x)
      f = -x(1) + 100*norm(x(2:end));
    end
    
    function f=fdiffpow(x)
      N = size(x,1); 
      if N < 2 
        error('dimension must be greater one'); 
      end
      f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));
    end
    
    function f=frastrigin10(x)
      N = size(x,1); 
      if N < 2 
        error('dimension must be greater one'); 
      end
      scale=10.^((0:N-1)'/(N-1));
      f = 10*size(x,1) + sum((scale.*x).^2 - 10*cos(2*pi*(scale.*x)));
    end
    
    function f=frand(~)
      f=rand;
    end
  end
end