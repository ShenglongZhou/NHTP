To implement this solver, please 

[1] run startup.m first to add the path;
[2] run demonXXXX.m to solve different problems

This source code contains the algorithm described in:
"S. Zhou, N. Xiu and H. Qi, Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit, 
Journal of Machine Learning Research, 22(12):1-45, 2021"

Please give credits to this paper if you use the code for your research.

% ===================================================================
% The citation of the solver  NHTP takes the form of
%               
%                 Out = NHTP(problem,data,n,s,pars)
% 
% It  aims at solving the sparsity constrained optimization with form
%
%         min_{x\in R^n} f(x)  s.t.  \|x\|_0<=s
%
% where s is the given sparsity, which is << n.  
%
% Inputs:
%     problem:  A text string for different problems to be solved, (required)
%               = 'CS',  compressed sensing problems
%               = 'LCP', linear complementarity problems
%               = 'LR',  sparse logistic regression problems
%               = 'SCO', other sparsity constrained optimization problems
%     data    : A triple structure, (required)
%               data.A, the measurement matrix, or a function handle @(x)A(x);
%               data.At = data.A',or a function handle @(x)At(x);
%               data.b, the observation vector 
%     n       : Dimension of the solution x, (required)
%     s       : Sparsity level of x, an integer between 1 and n-1, (required)
%     func    : function handle, define the function value, gradient, Hessian of f(x)
%               it has the form: [out1,out2] = func(x,flag,T1,T2)               
%     pars:     Parameters are all OPTIONAL
%               pars.x0      --  Starting point of x,   pars.x0=zeros(n,1) (default)
%               pars.eta     --  A positive parameter,  a default one is given related to inputs  
%               pars.display --  =1. Display results for each iteration.(default)
%                                =0. Don't display results for each iteration.
%               pars.draw    --  A  graph will be drawn if pars.draw=1 (default) 
%                                No graph will be drawn if pars.draw=0
%               pars.maxit   --  Maximum number of iterations, (default,2000) 
%               pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%
% Outputs:
%     Out.sol:           The sparse solution x
%     Out.sparsity:      Sparsity level of Out.sol
%     Out.normgrad:      L2 norm of the gradient at Out.sol  
%     Out.error:         Error used to terminate this solver 
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%     Out.obj:           Objective function value at Out.sol 

% =================================================================
% Example I:  compressed sensing problem

n         = 2000; 
m         = ceil(0.25*n);
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n);
x(I(1:s)) = randn(s,1);
data.A    = randn(m,n)/sqrt(n);
data.b    = data.A*x ;
data.At   = data.A';
pars.eta = 1;
out       = NHTP('CS',data,n,s,pars);
ReoveryShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example II:  linear complementarity problem 

n         = 2000; 
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n); I = I(1:s);
x(I)      = rand(s,1);
A         = randn(n,ceil(n/4));
data.A    = A*A'/n;  Ax=data.A*x;
data.b    = abs(Ax); data.b(I)=-Ax(I); 
data.At   = data.A';
pars.eta  = 1;
out       = NHTP('LCP',data,n,s,pars);
ReoveryShow(out.sol,x,[900,500,500,250],1)

% =================================================================
% Example III:  Logistic regression problem

n         = 2000; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     
I         = randperm(n);
I         = I(1:s); 
data.A    = randn(m,n); 
data.At   = data.A'; 
q         = 1./(1+exp(-data.A(:,I)*randn(s,1)));
data.b    = zeros(m,1);
for i     = 1:m    
data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
end               

pars.eta  = 0.5;
out       = NHTP('LR',data,n,s,pars);
