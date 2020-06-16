clc; clear; close all;

ExNam   = 1; %= 1, 2, 3, 4 or 5                          
ExMat   = 2; %= 1 or 2

n       = 10000;  
m       = ceil(0.25*n); 
s       = ceil(0.05*n);
 
switch ExNam
    
    case 1 % demon compressed sensing problems
    MatType = {'GaussianMat','PartialDCTMat'}; 
    data    = compressed_sensing_data(MatType{ExMat}, m,n,s,0);
    funcname='compressed_sensing'; 
    
    case 2 % demon sparse logistic regression problems
    MatType = {'Indipendent','Correlated'};
    rho     = 0.5; % 0<= rho <=1 it is useful for 'Correlated' data
    data    = logistic_random_data(MatType{ExMat},m,n,s,rho);
    funcname='logistic_regression'; 
    
    case 3 % demon sparse linear complementarity problems
    MatType = {'z-mat','sdp','sdp-non'};
    data    = lcp_data(MatType{ExMat},n,s);
    funcname='slcp'; pars.eta=1;
    
    case 4 % demon a simple example 
    a       = 0.1*randn; b=0.1*rand(n,1);
    data    = @(x,fgh)simple_ex2_func(x,fgh, a, b);
    funcname='general_example'; 
    
    case 5 % demon a simple example 
    data    = @(x,fgh)simple_ex4_func(x,fgh);
    funcname='general_example';  n=2; s=1;
end

pars.draw = 1;
fun       = str2func(funcname);
func      = @(x,fgh,T1,T2)fun(x,fgh,T1,T2,data);
out       = NHTP(n,s,func,pars) 


 
