function data = compressed_sensing_data(problemname,m,n,s,nf )
% This file aims at generating data of 3 examples, Gaussian, Partial DCT, 
% and Toeplitz Correlation type measurement matrices.
% Inputs:
%       problemname -- can be 'Gaussian','PartialDCTMat' or 'ToeplitzCorMat'
%       m and n     -- dimensions of A
%       s           -- sparsity level of x_opt, an intger between 1 and n-1
%       nf          -- noise ratio
% Outputs:
%       data.A           --  m x n order measurement matrices,  (required)
%       data.At          --  the transpose of data.A, i.e., data.At=data.A', (required)
%       data.b           --  m x 1 order observation vector, (required)
%       data.x_opt       --  n x 1 order 'true' sparse solution, (optional) 
%
% which satisfies       b = A*x_opt + nf*noise 

% written by Shenglong Zhou, 13/10/2018

start = tic;
fprintf('Start to generate the compressed sensing data...\n'); 

A = zeros(m,n);
 
switch problemname
     
    case 'GaussianMat'        
        A      = randn(m,n);
        I0     = randperm(n); 
        I      = I0(1:s);

    case 'PartialDCTMat'
        r      = rand(m,1);
        column = 1:n;
        for i  = 1:m 
        A(i,:) = cos(2*pi*r(i)*(column-1));
        end
        I0     = randperm(n); 
        I      = I0(1:s);

    case 'ToeplitzCorMat'
        Sig    = zeros(n,n);
        t      = 1:n; 
        for i  = 1:n
        Sig(i,:)= (.5).^(abs(i-t)); 
        end
        Sig    = real(Sig^(1/2));   
        A      = randn(m,n)*Sig;
        I0     = randperm(n);  
        I      = I0(1:s);   
        
    otherwise
        disp('input a problen name');        
end

x_opt     = zeros(n,1);  
while nnz(x_opt)~=s  
x_opt(I)  = randn(s,1); 
end
 
data.A     = normalization(A, 3);                   % required
data.b     = data.A(:,I)*x_opt(I)+nf*randn(m,1);    % required
data.At    = data.A';                               % required

data.x_opt = x_opt;                                 % optional

fprintf(' Data generation used %2.4f seconds.\n\n',toc(start)); 

end

