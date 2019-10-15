% demon compressed sensing problems with randomly generated data
clc; clear; close all;

n       = 256;  
m       = ceil(n/4); 
s       = ceil(0.01*n);                      
 
% % You could input any data including (data.A, data.At, data.b) 
data.A    = randn(m,n)/sqrt(m);
data.At   = data.A';
data.b    = randn(m,1)/sqrt(m);  

% Or you could input data from our data generation function
% ExMat = 1;
% MatType = {'GaussianMat','PartialDCTMat'}; 
% data    = compressed_sensing_data(MatType{ExMat}, m,n,s,0);

pars.eta  = 1;
func      = @(x,fgh,T1,T2)compressed_sensing(x,fgh,T1,T2,data);
out       = NHTP(n,s,func,pars) 

fprintf('\n Sample size:       m=%d,n=%d\n', m,n);
fprintf(' Recovery time:    %6.3fsec\n',  out.time);
if isfield(data,'x_opt')
fprintf(' Recovery accuracy: %5.2e\n\n', norm(out.sol-data.x_opt)/norm(data.x_opt));
end