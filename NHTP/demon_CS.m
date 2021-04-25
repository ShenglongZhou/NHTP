% demon compressed sensing problems 
clc; clear; close all;

n     = 10000;  
m     = ceil(n/4); 
s     = ceil(0.01*n);                      
test  = 1; 

switch test
  case 1       % Input data by our data generation function
       ExMat    = 1;
       MatType  = {'GaussianMat','PartialDCTMat'}; 
       nf       = 0;
       data     = compressed_sensing_data(MatType{ExMat},m,n,s,nf);
       if nf    > 0; pars.eta = 0.5; end
  case 2       % Input any data including (data.A, data.At, data.b), e.g.,
       I        = randperm(n); 
       Tx       = I(1:s);
       xopt     = zeros(n,1);  
       xopt(Tx) = randn(s,1); 
       data.A   = randn(m,n)/sqrt(m);
       data.At  = data.A';
       data.b   = data.A*xopt;  
       data.xopt= xopt; 
end
 
pars.draw = 0;
out       = NHTP('CS',data,n,s,pars);  
fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Objective:         %5.2e\n',  out.obj);
fprintf(' Sample size:       %dx%d\n', m,n);
if isfield(data,'xopt') && s<=1000
   ReoveryShow(data.xopt,out.sol,[1000, 550, 400 200],1)
end
