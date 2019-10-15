% demon sparse logistic regression problems with real data
clc; close all; clear;

prob      = 'newsgroup'; %'colon-cancer'
measure   = load(strcat(prob,'.mat')); 
label     = load(strcat(prob,'_label.mat'));   
label.b(label.b==-1)= 0;
[m,n]     = size(measure.A);  
data.A    = normalization(measure.A,1+(m>= 1000)); 
data.At   = data.A';
data.b    = label.b; 
clear measure label;

s         = ceil(0.01*m);
pars.eta  = 1;
func      = @(x,fgh,T1,T2)logistic_regression(x,fgh,T1,T2,data);
out       = NHTP(n,s,func);

saveas(figure(1), [pwd strcat(strcat('/outputs/',prob))]);  
saveas(figure(1), [pwd strcat(strcat('/outputs/',prob),'.eps')],'epsc');

