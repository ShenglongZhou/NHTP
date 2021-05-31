function [out1,out2] = CS(x,T1,T2,data)

if ~isa(data.A, 'function_handle') % A is a matrix 
    if  isempty(T1) && isempty(T2) 
        if  nnz(x) >= 0.8*length(x)
            Axb     = data.A*x-data.b;
        else
            Tx      = find(x); 
            Axb     = data.A(:,Tx)*x(Tx)-data.b;
        end
            out1    = (Axb'*Axb)/2;               % objective function value of f
        if  nargout == 2
            out2    = (Axb'*data.A)';                % gradien of f
        end
    else        
            AT = data.A(:,T1); 
        if  length(T1)<3000
            out1 = AT'*AT;                        %subHessian containing T1 rows and T1 columns
        else
            out1 = @(v)( (AT*v)'*AT )';      
        end       
        if  nargout == 2
            out2 = @(v)( (data.A(:,T2)*v)'*AT )'; %subHessian containing T1 rows and T2 columns
        end       
    end
else  % A is a function handle A*x=A(x)  
    if ~isfield(data,'At') 
        disp('The transpose-data.At-is missing'); return; 
    end
    if ~isfield(data,'n')  
        disp('The dimension-data.n-is missing');  return;  
    end   
    if  isempty(T1) && isempty(T2)  
        Axb  = data.A(x)-data.b;
        out1 = (Axb'*Axb)/2;              % objective function value of f
        if  nargout>1 
            out2 = data.At(Axb);          % gradien of f
        end
    else
        func = fgH(data);    
        out1 = @(v)func(v,T1,T1);         % subHessian containing T1 rows and T1 columns
        if  nargout>1
            out2 = @(v)func(v,T1,T2);     % subHessian containing T1 rows and T1 columns
        end  
        
    end
end

end

function Hess = fgH(data)
    suppz     = @(z,t)supp(data.n,z,t);
    sub       = @(z,t)z(t,:);
    Hess      = @(z,t1,t2)(sub( data.At( data.A(suppz(z,t2))),t1)); 
end

function z = supp(n,x,T)
    z      = zeros(n,1);
    z(T)   = x;
end



