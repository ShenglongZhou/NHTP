function demon_sco()
% demon a general sparsity constrained problem

clc; clear; close all;
pars.eta  = 0.1;
n         = 2;
s         = 1;
data      = @(x,fgh)simple_sco(x,fgh);  
out       = NHTP('SCO',data,n,s,pars);

%-----------------------------------------------------
function data = simple_sco(x,fgH)
    % This code provides information for
    %     min   x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1) 
    %     s.t. \|x\|_0<=s
    % where s=1
    a = sqrt(sum(x.*x)+1);
    switch fgH    
        case 'obj'  % objective function
        data = x'*[6 5;5 8]*x+[1 9]*x-a;        
        case 'grad'  % gradient 
        data = 2*[6 5;5 8]*x+[1; 9]-x./a;        
        case 'hess'  % Hessian matrix 
        data = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3;  
    end   
end

end
