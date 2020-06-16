% demon sparse logistic regression problems with real data
clc; close all; clear;

prob      = 'newsgroup'; %'colon-cancer'
Mat       = load(strcat(prob,'.mat')); 
label     = load(strcat(prob,'_label.mat'));   
label.b(label.b==-1)= 0;
[m,n]     = size(Mat.A);  
data.A    = normalization(Mat.A,1+(m>= 1000)); 
data.At   = data.A';
data.b    = label.b; 
clear Mat label;

s         = ceil(0.2*m);
pars.eta  = 5;
fname     = str2func('logistic_regression');
func      = @(x,fgh,T1,T2)fname(x,fgh,T1,T2,data); clear data
out       = NHTP(n,s,func) 


fprintf('\n Sample size:   m=%4d,n=%4d\n', m,n);
fprintf(' CPU time:      %.3fsec\n',  out.time);
fprintf(' Logistic Loss: %5.2e\n\n', out.obj);
path      = strcat(strcat('/outputs/',prob));
saveas(figure(1), [pwd path]);  
saveas(figure(1), [pwd path,'.eps'],'epsc');
