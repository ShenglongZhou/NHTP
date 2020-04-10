% demon compressed sensing problems with randomly generated data
clc; clear; close all;

n     = 10000;  
m     = ceil(n/4); 
s     = ceil(0.01*n);                      
test  = 2; 

switch test
    case 1   % Input any data including (data.A, data.At, data.b) 
    data.A   = randn(m,n)/sqrt(m);
    data.At  = data.A';
    data.b   = randn(m,1)/sqrt(m);  
    case 2   % Input data from our data generation function
	ExMat    = 1;
	MatType  = {'GaussianMat','PartialDCTMat'}; 
	data     = compressed_sensing_data(MatType{ExMat}, m,n,s,0);
    xopt     = data.x_opt;
    case 3   % Input real data
    load 'DrivFace.mat'; load 'nlab.mat'; %'identity.mat';
    [m,n]    = size(A);
    s        = ceil(0.2*min(m,n));
    data.A   = A/sqrt(m); clear A
    data.At  = data.A';
    data.b   = y/sqrt(m); clear y
end
 
pars.eta = 1;
fname    = str2func('compressed_sensing');
func     = @(x,fgh,T1,T2)fname(x,fgh,T1,T2,data); clear data
out      = NHTP(n,s,func,pars);  

fprintf('\n Sample size:       %dx%d\n', m,n);
fprintf(' Recovery time:    %6.3fsec\n',  out.time);
fprintf(' Objective:         %5.2e\n',  out.obj);
if exist('xopt')
fprintf(' Recovery accuracy: %5.2e\n\n', ...
          norm(out.sol-xopt)/norm(xopt));
end
