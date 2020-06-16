% demon sparse logistic regression problems with randomly generated data
clc; clear; close all;

n    = 10000;  
m    = ceil(n/5); 
s    = ceil(0.05*n);
test = 2;                     
switch test
    case 1  % You could input any data including (data.A, data.At, data.b), e.g, 
         I0      = randperm(m);  
         I       = I0(1:ceil(m/2)); 
         b       = ones(m,1);    
         b(I)    = 0;
         data.A  = repmat(b.*rand(m,1),1,n) + randn(m,n);
         data.At = data.A';    
         data.b  = b;
    case 2  % Or you could input data by our data generation function 
         ExMat   = 2;
         MatType = {'Correlated','Weakly-Indipendent'};
         data    = logistic_random_data(MatType{ExMat},m,n,s,0.5); 
end

pars.eta = 5;
fname    = str2func('logistic_regression');
func     = @(x,fgh,T1,T2)fname(x,fgh,T1,T2,data); clear data
out      = NHTP(n,s,func); 

fprintf('\n Sample size:   m=%4d,n=%4d\n', m,n);
fprintf(' CPU time:      %.3fsec\n',  out.time);
fprintf(' Logistic Loss: %5.2e\n\n', out.obj);
