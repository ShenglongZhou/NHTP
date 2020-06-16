function [out1,out2] = compressed_sensing(x,fgh,T1,T2,data)
    Tx = find(x); 
    Ax = data.A(:,Tx)*x(Tx);
    switch fgh        
    case 'ObjGrad'
        Axb  = Ax-data.b;
        out1 = sum(Axb.*Axb)/2;                %objective function value of f
        if  nargout>1 
        out2 = data.At*Axb;                    %gradien of f
        end
    case 'Hess'
        
        if length(T1)<5e3
        out1 = data.At(T1,:)*data.A(:,T1);     %submatrix containing T1 rows and T1 columns of Hessian
        else
        out1 = @(v)(data.At(T1,:)*(data.A(:,T1)*v));      
        end
        
        if nargout>1
        out2 = @(v)(data.At(T1,:)*(data.A(:,T2)*v));     %submatrix containing T1 rows and T2 columns of Hessian
        end  
        
    end
end


