function [out1,out2] = sco(x,fgh,T1,T2,data)
    switch fgh
    case 'ObjGrad'             
        out1 = data(x,'obj');  % objective function   
        if nargout>1
        out2 = data(x,'grad');  % gradient
        end       
    case 'Hess'
        H    = data(x,'hess'); 
        out1 = H(T1,T1);       % submatrix containing T1 rows and T1 columns of Hessian
        if nargout>1 
        out2 = H(T1,T2);       % submatrix containing T1 rows and T2 columns of Hessian
        end
        clear  H;
    end     
end



