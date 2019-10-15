function out = HTPCP(n, func, pars)
% This code aims at solving the sparse complementarity problem (CL) with form
%
%                  min    sum_i |x_i|^0.5
%                  s.t.   x = max(0, x-F(x))
%
% where x\in\R^n, F(x):\R^n-->\R^n.  
%
% Inputs:
%     n       : Dimension of the solution x, (required)
%     func    : function handle defines the function F(x). For example,
%               A linear CP:    F(x)=M*x+q, where M\R^{n*n},q\R^{n}
%                         func = @(x,T)(M(:,T)*x(T)+q);                
%               A nonlinear CP: F(x)=a.*arctan(x)+M*x+q, where a\R^{n},M\R^{n*n},q\R^{n}
%                         func = @(x,T)(a.*atan(x)+M(:,T)*x(T)+q); 
%               
%     pars:     Parameters are all OPTIONAL
%               pars.iteron --  Results will  be shown for each iteration if pars.iteron=1 (default)
%                               Results won't be shown for each iteration if pars.iteron=0 
%               pars.maxit  --  Maximum nonumber of iteration.  pars.maxit=5000 (default) 
%               pars.tol    --  Tolerance of stopping criteria. pars.maxit=1e-5 (default) 
%
% Outputs:
%     Out.x:             The sparse solution x 
%     Out.Fx:            F(x)
%     Out.xFx:           x'*F(x)
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%
%  This solver was created based on the algorithm proposed by  
%  Shang, M.,Zhang, C.,Peng, D.,& Zhou, S., A half thresholding projection 
%  algorithm for sparse solutions of LCPs, Optimization Letters,9(6),1231-1245,2015
%
%%%%%%%    Send your comments and suggestions to                     %%%%%%
%
%%%%%%%    shenglong.zhou@soton.ac.uk                                %%%%%%
% 
%%%%%%%    Warning: Accuracy may not be guaranteed!!!!!              %%%%%%
warning off;

if nargin<2; error('Imputs are not enough!\n'); end
if nargin<3; pars=[]; end
if isfield(pars,'iteron');iteron = pars.iteron; else; iteron = 1;     end
if isfield(pars,'maxit'); maxit  = pars.maxit;  else; maxit  = 5000;  end
if isfield(pars,'tol');   tol    = pars.tol;    else; tol    = 1e-5;  end  
   
t0       = tic;
tau      = .75;
lambda   = 10; 
t        = 2/3; 
z        = zeros(n,1);
x        = zeros(n,1); 

if iteron 
fprintf(' Start to run the sover...\n'); 
fprintf('\n Iter          Error          Time \n'); 
fprintf('------------------------------------\n');
end

for iter = 1: maxit
    
    zxk  = norm(x-z);
    
    if iter>1
       err = zxk /max(1,norm(x(T)));
       if mod(iter,5)==0 && iteron
       fprintf('%4d         %5.2e       %5.2fsec\n',iter,err,toc(t0)); 
       end
       if err<tol; break; end
    end  
    
    %--------------------------- x_k1 iteration---------------------------%            
    if mod(iter,10)==0
       lambda = max(1e-10,tau^(iter/100)*lambda); 
    end
     x0     = x;
     x      = zeros(n,1);
     T      = find( z > 54^(1/3)/4*lambda^t );
     x(T)   = t*z(T).*( 1+cos( t*pi-t*acos( lambda./( abs(z(T))/3 ).^t/8 ) ) );
     %--------------------------- z_k iteration---------------------------%
     alpha = 1;
     f1    = norm(x-x0)^2 + zxk^2;
     f2    = norm(x-z)^2;
     Fx    = func(x,T);                  
     for j = 1:5
     z     = max(0,x-alpha*Fx);
     if norm(x-z)^2 + alpha*f1 < f2; break; end
     alpha = alpha/2;     
     end

end

fprintf('------------------------------------\n');
out.x    = x;
out.Fx   = Fx;
out.xFx  = sum(x(T).*Fx(T));
out.iter = iter;
out.time = toc(t0);
end
