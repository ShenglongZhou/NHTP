To implement this solver, please 

[1] run startup.m first to add the path;
[2] run demonXXXX.m to solve different problemsa

This source code contains the algorithm described in:
"S. Zhou, N. Xiu and H. Qi, Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit, 
Journal of Machine Learning Research, 22(12):1-45, 2021"

Please give credits to this paper if you use the code for your research.

% ===================================================================
% The citation of the solver  NHTP takes the form of
%               
%                 Out = NHTP(func,n,s,pars)
% 
% It  aims at solving the sparsity constrained optimization with form
%
%         min_{x\in R^n} f(x)  s.t.  \|x\|_0<=s
%
% where s is the given sparsity, which is << n.  
%
% Inputs:
%     func:   A function handle defines (objective,gradient,sub-Hessain)  (required)
%     n       : Dimension of the solution x,                              (required)
%     s       : Sparsity level of x, an integer between 1 and n-1,        (required)    
%           
%     pars:     Parameters are all OPTIONAL
%               pars.x0      --  Starting point of x,   pars.x0=zeros(n,1) (default)
%               pars.eta     --  A positive parameter,  a default one is given related to inputs  
%               pars.display --  Display results or not for each iteration (default, 1)
%               pars.draw    --  Draw or not draw a graph (default, 0) 
%               pars.maxit   --  Maximum number of iterations, (default,2000) 
%               pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%
% Outputs:
%     Out.sol:           The sparse solution x
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%     Out.obj:           Objective function value at Out.sol 


% Here are some examples that you can run
% =================================================================
% Example I:  compressed sensing problem

n         = 10000; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     
x         = zeros(n,1);
I         = randperm(n); 
I         = I(1:s);
x(I)      = randn(s,1);
data.A    = randn(m,n)/sqrt(m);
data.b    = data.A(:,I)*x(I);
func      = @(x,T1,T2)CS(x,T1,T2,data);
out       = NHTP(func,n,s); 
ReoveryShow(x,out.sol,[900,500,500,250],1)

% =================================================================
% Example II:  linear complementarity problem 

n         = 10000; 
s         = ceil(0.01*n);     
x         = zeros(n,1);
I         = randperm(n); 
T         = I(1:s);
x(T)      = rand(s,1);
A         = randn(n,ceil(n/4));
data.A    = A*A'/n;  
data.At   = data.A';
Ax        = data.A(:,T)*x(T);
data.b    = abs(Ax); 
data.b(T) = -Ax(T); 
func      = @(x,T1,T2)LCP(x,T1,T2,data);
out       = NHTP(func,n,s);
ReoveryShow(x,out.sol,[900,500,500,250],1)

% =================================================================
% Example III:  Logistic regression problem

n         = 10000; 
m         = ceil(0.25*n);
s         = ceil(0.05*n);     
I         = randperm(n);
T         = I(1:s); 
data.A    = randn(m,n); 
q         = 1./(1+exp(-data.A(:,T)*randn(s,1)));
data.b    = zeros(m,1);
for i     = 1:m    
data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
end               
func      = @(x,T1,T2)LogitReg(x,T1,T2,data);
out        = NHTP(func,n,s);
