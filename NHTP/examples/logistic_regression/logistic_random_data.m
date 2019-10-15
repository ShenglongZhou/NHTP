function [data,out] = logistic_random_data(type,m,n,s,rho)

% This file aims at generating data of 3 examples
% Inputs:
%       problemname -- can be 'Indipendent','Correlated' or 'Weakly-Indipendent' 
%       m and n     -- dimensions of A
%       s           -- sparsity level of x, an intger between 1 and n-1
%       nf          -- noise ratio
%       rho         -- a scalar from [0,1]. This is typical for 'Corrolated' measurement matrix
% Outputs:
%       data.A      --  m x n order measurement matrices, (required)
%       data.At     --  the transpose of data.A, i.e., data.At=data.A', (required)
%       data.b      --  m x 1 order response vector, (required) 
%       data.x      --  n x 1 order vector used for generating the response y, (optional)
%
% written by Shenglong Zhou, 13/10/2018

start = tic;
fprintf('Start to generate the logistic regression data...\n'); 
switch type
    case 'Indipendent'
        I0    = randperm(m);
        b     = ones(m,1);
        I     = I0(1:ceil(m/2)); 
        b(I)  = 0;
        A     = repmat(b.*rand(m,1),1,n)+ randn(m,n);
        x     = [];
        out   = [];
    case 'Correlated'
        I0    = randperm(n);
        x     = zeros(n,1);
        I     = I0(1:s); 
        x(I)  = randn(s,1);

        v     = randn(m,n);
        A     = zeros(m,n); 
        A(:,1)= randn(m,1);

        for j=1:n-1
        A(:,j+1)=rho*A(:,j)+sqrt(1-rho^2)*v(:,j);
        end
        Ax    = A(:,I)*x(I);
        q     = 1./(1+exp(-Ax));
        
        b     = zeros(m,1);
        for i = 1:m    
        b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
        end
             
        out.ser = sum(abs(b-max(0,sign(Ax))))/m;
        out.f  = sum( log(1+exp(Ax))- b.*Ax )/m;
        
    case 'Weakly-Indipendent'       
        I0    = randperm(n);
        x     = zeros(n,1);
        I     = I0(1:s); 
        x(I)  = randn(s,1);
        
        A     = randn(m,n); 
        Ax    = A(:,I)*x(I);
        q     = 1./(1+exp(-Ax));
        b     = zeros(m,1);
        for i = 1:m    
        b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
        end
 
        out.ser = sum(abs(b-max(0,sign(Ax))))/m;
        out.f  = sum( log(1+exp(Ax))- b.*Ax )/m;
end  

data.A  = A;                   % required
data.b  = b;                   % required
data.At = A';                  % required

data.x  = x;                   % optional

fprintf(' Data generation used %2.4f seconds.\n\n',toc(start)); 

end
