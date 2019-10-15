% demon sparse complementarity problems by two phases procedure
clc; clear; close all;                    
n       = 2000;  
s       = ceil(0.01*n);
ExMat   = 2; %= 1, 2, 3

MatType  = {'z-mat','sdp','sdp-non'};
data     = lcp_data(MatType{ExMat},n,s);

% First phase
pars0.tol  = 1e-3;
fun        = @(x,T)(data.M(:,T)*x(T)+data.q);
out0       = HTPCP(n, fun , pars0)  

% Second phase
pars.x0  = out0.x;
func     = @(x,fgh,T1,T2)slcp(x,fgh,T1,T2,data);
out      = NHTP(n,s,func,pars)  

 
