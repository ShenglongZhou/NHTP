function data = lcp_data(example,n,s )
% This file aims at generating data of 3 examples, Gaussian, Partial DCT, 
% and Toeplitz Correlation type measurement matrices.
% Inputs:
%       example     -- can be 'z-mat','sdp' or 'sdp-non'
%       n           -- dimension  of M\in \R^{n x n}
%       s           -- sparsity level of xopt, an intger between 1 and n-1
% Outputs:
%       data.M      --  n x n order measurement matrices,  (required)
%       data.Mt     --  the transpose of data.M, i.e., data.Mt=data.M', (required)
%       data.q      --  n x 1 order observation vector, (required)
%       data.n      -- dimension  of M\in \R^{n x n}
%       data.xopt   --  n x 1 order 'true' sparse solution, (optional) 
%
% which satisfies       b = A*x_opt + nf*noise 

% written by Shenglong Zhou, 05/12/2019

start = tic;
fprintf(' Start to generate the linear complementarity problem data...\n'); 
switch example
    case 'z-mat'
         M       = eye(n)-ones(n)/n;
         q       = ones(n,1)/n; 
         q(1)    = 1/n-1;
         xopt    = zeros(n,1); 
         xopt(1) = 1;
         Mt      = M;
    case 'sdp'
         Z       = randn(n,ceil(n/2));
         M       = Z*Z';
         [xopt,T]= get_sparse_x(n,s); 
         Mx      = M(:,T)*xopt(T);
         q       = abs(Mx);
         q(T)    = -Mx(T); 
         Mt      = M;
         M=M/n;  Mt=Mt/n; q=q/n;
    case 'sdp-non'
         Z       = rand(n,ceil(n/4));
         M       = Z*Z';
         [xopt,T]= get_sparse_x(n,s); 
         Mx      = M(:,T)*xopt(T);
         q       = rand(n,1);
         q(T)    = -Mx(T); 
         Mt      = M;
         M=M/n;  Mt=Mt/n; q=q/n;      
end
    
data.M    = M;
data.Mt   = Mt;
data.q    = q;
data.n    = n;
data.xopt = xopt;
clear M Mt xopt Mx q Z
fprintf(' Data generation used %2.4f seconds.\n\n',toc(start)); 
end


function [x,T]=get_sparse_x(n,s)
I     = randperm(n); 
T     = I(1:s); 
x     = zeros(n,1); 
x(T)  = 0.1+abs(randn(s,1)); 
end
