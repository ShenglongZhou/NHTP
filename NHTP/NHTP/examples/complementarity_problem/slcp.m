function [out1,out2] = slcp(x,fgh,T1,T2,data)
    
    if isfield(data,'r') 
    r   = data.r ;
    else
    r   = 2;
    end
    n   = data.n;
    M   = data.M;
    Mt  = data.Mt;
    q   = data.q;     
    clear data; 
    
    eps =  0;
    ip  = find(x>eps);
    in  = find(x<-eps);
    ix  = union(ip,in); 
    Mx  = M(:,ix)*x(ix) + q;    
    tp  = find(Mx>eps);  
    tn  = find(Mx<-eps);  
    tt  = intersect(ip,tp); 
 
    if isempty(tt); com=0;else; com=1; end
    Mxn = abs(Mx(tn));     
    xn  = abs(x(in)); 
    
    switch fgh        
    case 'ObjGrad'
         
        out1 =  (sum(xn.^r) + sum(Mxn.^r) )/r ;                 %objective function    
        if com
        out1 =  out1 + sum( ( x(tt).* Mx(tt) ).^r )/r;  
        end
        
        if  nargout>1 
        out2     = - Mt(:,tn)* ( Mxn.^(r-1) );
        out2(in) = out2(in) - xn.^(r-1) ; 
        
        if com
        out2     = out2 +  Mt(:,tt)*( x(tt).^r.*(Mx(tt).^(r-1)) ); 
        out2(tt) = out2(tt) + ( x(tt).^(r-1) ).*( Mx(tt).^r ); 
        end
        
        end
        
    case 'Hess'
        s1   = length(T1);
        mx   = max(x,0);       
        mMx  = max(Mx,0);
        
        if r~=2
            r1   = r-1; r2 = r-2;          
            z2   = r1* ( Mxn.^r2 ) ;                
            xy   = r1*( ( mx(T1).^r2 ).*( mMx(T1).^r ) + (abs(min(x(T1),0))).^r2 );
            MM   = Mt(T1,tn)*( repmat(z2,1,s1).*M(tn,T1) );
            if com
            z1   = r1* ( mx(tt).^r.*( mMx(tt).^r2 ) );  
            MM   = MM+ Mt(T1,tt)*( repmat(z1,1,s1).*M(tt,T1) );    
            end
        else  
        	            
            z    = ones(n,1); z(ip) = mMx(ip).^2;
            xy   = z(T1);   
            tn0  = setdiff(1:n,tp);
            MM   = Mt(T1,tn0)*M(tn0,T1) ;
            if com
            z1   = mx(tt).^r ;
            MM   = MM+Mt(T1,tt)*(repmat(z1,1,s1).*M(tt,T1) );    
            end
        end  
        
        tem1 = r*( mx(T1).*mMx(T1) ).^(r-1);     
        %submatrix  containing T1 rows and T1 columns of Hessian
        out1 = repmat(tem1,1,s1).*M(T1,T1);
        out1 = out1 + out1' + MM;         
        out1(1:(s1+1):end) = out1(1:(s1+1):end) + xy';    
        
        if nargout>1 
        s2   = length(T2) ;
        if r~=2
        MM   = Mt(T1,tn)*(repmat(z2,1,s2).*M(tn,T2));
        else
        MM   = Mt(T1,tn0)*M(tn0,T2);    
        end
        
        if com
        MM   = MM + Mt(T1,tt)*(repmat(z1,1,s2).*M(tt,T2));   
        end
        
        tem2 = r*( mx(T2).*mMx(T2) ).^(r-1);     
         %submatrix  containing T1 rows and T2 columns of Hessian
        out2 = repmat(tem1,1,s2).*M(T1,T2) + Mt(T1,T2).*repmat(tem2',s1,1) + MM; 
        end
  
    end

    clear ip in tp tn;
end


