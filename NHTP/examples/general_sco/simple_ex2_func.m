function data = simple_ex2_func(x,fgH,a,b)

% This code provides information for
%     min  (1/n) sum_{i=1}^n { [(x_i-b_i)^2-a]^2/4 + sin(x_i) }, 
%     s.t. \|x\|_0<=s
 
n  = length(x);
switch fgH    
    case 'obj'  % objective function
    data = (sum((x-b).^2)-a)^2/4+sum(sin(x));     
    
    case 'grad'  % gradient 
    data = (sum((x-b).^2)-a)*(x-b)+cos(x);    
    
    case 'hess'  % Hessian matrix 
    data = 2*((x-b)*(x-b)')+sum((x-b).^2)*eye(n)-diag(cos(x));  
end

data = data /n;
end


 