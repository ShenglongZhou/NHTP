% demon a general sparsity constrained problem
%     min    x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1)  
%     s.t. \|x\|_0<=s
% where s=1
% you can find this function in 'examples'-->'general_sco'

clc; clear; close all;
pars.eta  = 0.1;
data      = @(x,fgh)simple_ex(x,fgh);  
out1      = NHTP('SCO',data,2,1,pars);
