function [out1,out2] = compressed_sensing(x,fgh,T1,T2,data)

if ~isa(data.A, 'function_handle')             % A is a matrix 
    Tx = find(x); 
    Ax = data.A(:,Tx)*x(Tx);
    switch fgh        
    case 'ObjGrad'
        Axb  = Ax-data.b;
        out1 = sum(Axb.*Axb)/2;                % objective function value of f
        if  nargout>1 
        out2 = data.At*Axb;                    % gradien of f
        end
    case 'Hess'
        
        if  length(T1)<2000
            out1 = data.At(T1,:)*data.A(:,T1);     %subHessian containing T1 rows and T1 columns
        else
            AT1  = data.A(:,T1);
            out1 = @(v)((AT1*v)'*AT1)';      
        end
        
        if nargout>1
        out2 = @(v)(data.At(T1,:)*(data.A(:,T2)*v)); %subHessian containing T1 rows and T2 columns
        end  
        
    end
else                                       % A is a function handle A*x=A(x)
    func = fgH(data);
    switch  fgh        
      case 'ObjGrad'
            out1 = func.obj(x);                % objective function value of f
            if  nargout>1 
            out2 = func.grad(x);               % gradien of f
            end
      case 'Hess'
            out1 = @(v)func.Hess(v,T1,T1);     % subHessian containing T1 rows and T1 columns
            if nargout>1
            out2 = @(v)func.Hess(v,T1,T2);     % subHessian containing T1 rows and T1 columns
            end  
        
    end
end

end

function func = fgH(data)
Axb       = @(z)data.A(z)-data.b;
func.obj  = @(z)norm(Axb(z))^2/2;              
func.grad = @(z)data.At(Axb(z));  
suppz     = @(z,t)supp(data.n,z,t);
sub       = @(z,t)z(t,:);
func.Hess = @(z,t1,t2)(sub( data.At( data.A(suppz(z,t2))),t1)); 
end

function z=supp(n,x,T)
z    = zeros(n,1);
z(T) = x;
end



