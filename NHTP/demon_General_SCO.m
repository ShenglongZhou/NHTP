% demon a general sparsity constrained problem
%     min    x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1)  
%     s.t. \|x\|_0<=s
% where s=1

clc; clear; close all;
n    = 2;    
s    = 1;

% you can find this function in 'examples'-->'general_sco'
pars.eta  = .1;
data      = @(x,fgh)simple_ex4_func(x,fgh);  
func      = @(x,fgh,T1,T2)general_example(x,fgh,T1,T2,data);
out1      = NHTP(n,s,func,pars)  

% or you can solve by following way, 
% where @exfunc can be found in the folder 'general_sco'
pars.eta  = .1;
out2      = NHTP(n,s,@exfunc,pars) 

if isfield(pars,'draw') && pars.draw
saveas(figure(1), 'outputs\GenSCO.eps','epsc');
saveas(figure(1), 'outputs\GenSCO.fig');
end
