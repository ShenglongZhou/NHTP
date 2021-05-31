% demon sparse logistic regression problems 
clc; clear; close all;

test = 1;
switch test
    case 1  % Or you could input data by our data generation function 
         n         = 10000;  
         m         = ceil(n/5); 
         s         = ceil(0.05*n);
         rho       = 0.5;
         I         = randperm(n);
         T         = I(1:s); 
         data.A    = randn(m,n); 
         q         = 1./(1+exp(-data.A(:,T)*randn(s,1)));
         data.b    = zeros(m,1);
         for i     = 1:m    
         data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]);
         end 
    case 2  % Or you could input real data including (data.A, data.b), e.g,
         prob      = 'colon-cancer'; 
         Mat       = load(strcat(prob,'.mat')); 
         label     = load(strcat(prob,'_label.mat'));   
         label.b(label.b==-1)= 0;
         [m,n]     = size(Mat.A);  
         s         = ceil(0.2*m);
         data.A    = normalization(Mat.A,1+(m>= 1000)); 
         data.b    = label.b; 
end
        
func = @(x,T1,T2)LogitReg(x,T1,T2,data);
out  = NHTP(func,n,s);
fprintf(' CPU time:       %.3fsec\n',  out.time);
fprintf(' Logistic Loss:  %5.2e\n', out.obj);
fprintf(' Sample size:    %dx%d\n', m,n);
