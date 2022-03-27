function data = CSdata(problemname,m,n,s,nf )
% This file aims at generating data of 3 examples, Gaussian, Partial DCT, 
% and Toeplitz Correlation type measurement matrices.
% Inputs:
%       problemname -- can be 'Gaussian','PartialDCTMat' or 'ToeplitzCorMat'
%       m and n     -- dimensions of A
%       s           -- sparsity level of xopt, an intger between 1 and n-1
%       nf          -- noise ratio
% Outputs:
%       data.A           --  m x n order measurement matrices,  (required)
%       data.b           --  m x 1 order observation vector, (required)
%       data.xopt        --  n x 1 order 'true' sparse solution, (optional) 
%
% which satisfies       b = A*xopt + nf*noise 

% written by Shenglong Zhou, 13/10/2018

start = tic;
fprintf(' Please wait for CS data generation ...\n'); 

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

xopt     = zeros(n,1);
while nnz(xopt)~= s  
      xopt(I)   = randn(s,1); 
end
xopt(I)    = xopt(I) + 2*nf*sign(xopt(I));
data.A     = normalization(A, 3);                   % required
data.b     = data.A(:,I)*xopt(I)+nf*randn(m,1);     % required
data.xopt  = xopt;                                  % optional
clear A xopt I0 I
fprintf(' Data generation used %2.4f seconds.\n\n',toc(start)); 

end
