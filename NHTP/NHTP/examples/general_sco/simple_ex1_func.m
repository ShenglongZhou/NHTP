function data = simple_ex1_func(x,fgH,a)

% This code provides information for
%     min  (1/2n) sum_{i=1}^n (x_i-a_i-x_1*...*x_n)^2 
%     s.t. \|x\|_0<=s
 
n  = length(x);

switch fgH
    
case 'obj'   % objective function
      data  = 0.5*sum((x-a-prod(x)).^2);  
      
case 'grad'  % gradient  
      x(x==0)=1e-10;
      px   = prod(x);
      data  = (x-a-prod(x)).*(1+px./x);   
      
case 'hess' % Hessian matrix  
      x(x==0)=1e-10;
      px  = prod(x);
      data = (px./(x*x')).*(repmat((2*(x+px)-a),1,n));
      data(1:n+1:end)=(1+px./x).^2;       
end
 
data = data/n; 


end


 