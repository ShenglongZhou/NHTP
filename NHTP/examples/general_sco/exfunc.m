function [out1,out2] = exfunc(x,fgh,T1,T2)  
a = sqrt(sum(x.*x)+1); 
switch fgh 
    case 'ObjGrad'  
        out1 = x'*[6 5;5 8]*x+[1 9]*x-a;  
        if nargout>1 
        out2 = 2*[6 5;5 8]*x+[1; 9]-x./a; 
        end  
    case 'Hess'  
        H = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3;  
        out1 = H(T1,T1); 
        if nargout>1  
        out2 = H(T1,T2);  
        end  
        clear H;
end 
end 