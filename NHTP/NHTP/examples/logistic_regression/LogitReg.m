function [out1,out2] = LogitReg(x,T1,T2,data)
  
   if  nnz(x) >= 0.8*length(x)
       Ax  = data.A*x; 
   else
       Tx  = find(x); 
       Ax  = data.A(:,Tx)*x(Tx); 
   end
    m      = length(data.b);
    eAx    = exp(Ax);
    mu     = 1e-6/m; 
        
    if  isempty(T1) && isempty(T2) 
        
        if  sum(eAx)==Inf 
            Tpos = find(Ax>300); 
            Tneg = setdiff(1:m,Tpos);
            obj  = sum(log(1+eAx(Tneg)))+sum(Ax(Tpos))-sum(data.b.*Ax);                
        else
            obj  = sum(log(1+eAx)-data.b.*Ax); 
        end
        out1 = obj/m;                                     %objective function 
        
        if  nargout == 2 
            out2  = ((1-data.b-1./(1+eAx))'*data.A)'/m + mu*x;         %gradien
        end
        
    else 
        eXx  = 1./(1+eAx);
        d    = eXx.*(1-eXx)/m; 
        XT   = data.A(:,T1);
        s    = length(T1);
        if s < 2000
            out1 = (d.*XT)'*XT + mu*speye(s);           %submatrix  containing T1 rows and T1 columns of Hessian
        else
            out1 = @(v)( mu*v+( (d.*(XT*v))'*XT )' );
        end
        if nargout == 2
            out2 = @(v)( (d.*(data.A(:,T2)*v))'*XT )';  %submatrix  containing T1 rows and T2 columns of Hessian   
        end
    end
     
end



