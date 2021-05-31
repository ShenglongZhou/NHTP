% demon a general sparsity constrained problem
function demon_SCO()
clc; clear; close all;
pars.eta  = 0.1;
n         = 2;
s         = 1; 
out       = NHTP(@func,n,s,pars);
fprintf(' CPU time:       %.3fsec\n',  out.time);
fprintf(' Objective:      %.4f\n', out.obj); 

%---------------------------------------------------------
function [out1,out2] = func(x,T1,T2)
    % This code provides information for
    %     min   x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1) 
    %     s.t. \|x\|_0<=s
    % where s=1
    a   = sqrt(sum(x.*x)+1);
    if  isempty(T1) && isempty(T2)    
        out1 = x'*[6 5;5 8]*x+[1 9]*x-a;        
        if  nargout == 2
            out2 = 2*[6 5;5 8]*x+[1; 9]-x./a;  
        end
    else
        H    = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3;  
        out1 = H(T1,T1);
        if  nargout == 2
            out2 = H(T1,T2);  
        end
    end   
end

end
