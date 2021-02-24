% demon sparse complementarity problems by two phases procedure
clc; clear; close all;                    
n         = 10000;  
s         = ceil(0.01*n);
ExMat     = 2; %= 1, 2, 3

MatType   = {'z-mat','sdp','sdp-non'};
data      = lcp_data(MatType{ExMat},n,s);
pars.draw = 0; 
out       = NHTP('LCP',data,n,s,pars);

fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Objective:         %5.2e\n',  out.obj);
fprintf(' Sample size:       %dx%d\n', n,n);
if isfield(data,'xopt') && s<=100
   ReoveryShow(data.xopt,out.sol,[1000, 550, 400 200],1)
end
