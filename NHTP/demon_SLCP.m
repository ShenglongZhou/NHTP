% demon sparse complementarity problems by two phases procedure
clc; clear; close all;                    
n       = 10000;  
s       = ceil(0.01*n);
ExMat   = 2; %= 1, 2, 3

MatType  = {'z-mat','sdp','sdp-non'};
data     = lcp_data(MatType{ExMat},n,s);

% First phase (one can remove this phase)
pars0.tol  = 1e-3;
fun        = @(x,T)(data.M(:,T)*x(T)+data.q);
out0       = HTPCP(n, fun , pars0);  

% Second phase
if  exist('out0')
    pars.x0  = out0.x;
else
    pars.x0 = zeros(n,1);
end
func     = @(x,fgh,T1,T2)slcp(x,fgh,T1,T2,data);
out      = NHTP(n,s,func,pars);  

fprintf(' CPU time:          %.3fsec\n',  out.time);
fprintf(' Objective:         %5.2e\n',  out.obj);
fprintf(' Sample size:       %dx%d\n', n,n);
if isfield(data,'xopt') && s<=100
   ReoveryShow(data.xopt,out.sol,[1000, 550, 400 200],1),pause(0.05)
end
