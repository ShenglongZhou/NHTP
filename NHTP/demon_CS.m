% demon compressed sensing problems 
clc; clear; close all;

n     = 10000;  
m     = ceil(n/4); 
s     = ceil(0.01*n);                      
test  = 1; 

switch test
  case 1       % Input data by our data generation function
       ExMat   = 1;
       MatType = {'GaussianMat','PartialDCTMat'}; 
       nf      = 0;
       data    = compressed_sensing_data(MatType{ExMat},m,n,s,nf);
       pars.eta= 0.5*(nf>0)+(nf==0);
  case 2       % Input any data including (data.A, data.At, data.b), e.g.,
       data.A  = randn(m,n)/sqrt(m);
       data.At = data.A';
       data.b  = randn(m,1)/sqrt(m);  
       pars.eta= 1;
end
 
pars.draw = 0;
out       = NHTP('CS',data,n,s,pars);  
fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Objective:         %5.2e\n',  out.obj);
fprintf(' Sample size:       %dx%d\n', m,n);
if isfield(data,'xopt') && s<=100
   ReoveryShow(data.xopt,out.sol,[1000, 550, 400 200],1)
end
