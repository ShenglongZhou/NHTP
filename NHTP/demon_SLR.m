% demon sparse logistic regression problems 
clc; clear; close all;

test = 2;

switch test
    case 1  % Or you could input data by our data generation function 
         ExMat   = 2;
         MatType = {'Correlated','Weakly-Indipendent'};
         n       = 10000;  
         m       = ceil(n/5); 
         s       = ceil(0.05*n);
         data    = logistic_random_data(MatType{ExMat},m,n,s,0.5); 
    case 2  % Or you could input real data including (data.A, data.At, data.b), e.g,
         prob      = 'colon-cancer'; 
         Mat       = load(strcat(prob,'.mat')); 
         label     = load(strcat(prob,'_label.mat'));   
         label.b(label.b==-1)= 0;
         [m,n]     = size(Mat.A);  
         s         = ceil(0.2*m);
         data.A    = normalization(Mat.A,1+(m>= 1000)); 
         data.At   = data.A';
         data.b    = label.b; 
         clear Mat label;
end

out = NHTP('LR',data,n,s);
fprintf(' CPU time:       %.3fsec\n',  out.time);
fprintf(' Logistic Loss:  %5.2e\n', out.obj);
fprintf(' Sample size:    %dx%d\n', m,n);
