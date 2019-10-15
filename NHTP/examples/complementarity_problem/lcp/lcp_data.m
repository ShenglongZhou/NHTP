function data = lcp_data(example,n,s )
switch example
    case 'z-mat'
         M       = eye(n)-ones(n)/n;
         q       = ones(n,1)/n; 
         q(1)    = 1/n-1;
         xopt    = zeros(n,1); 
         xopt(1) = 1;
         Mt      = M;
    case 'sdp'
         Z       = randn(n,ceil(n/2));
         M       = Z*Z';
         [xopt,T]= get_sparse_x(n,s); 
         Mx      = M(:,T)*xopt(T);
         q       = abs(Mx);
         q(T)    = -Mx(T); 
         Mt      = M;
         M=M/n;  Mt=Mt/n; q=q/n;
    case 'sdp-non'
         Z       = rand(n,ceil(n/4));
         M       = Z*Z';
         [xopt,T]= get_sparse_x(n,s); 
         Mx      = M(:,T)*xopt(T);
         q       = rand(n,1);
         q(T)    = -Mx(T); 
         Mt      = M;
         M=M/n;  Mt=Mt/n; q=q/n;      
end
    
data.M    = M;
data.Mt   = Mt;
data.q    = q;
data.n    = n;
data.xopt = xopt;
end


function [x,T]=get_sparse_x(n,s)
I     = randperm(n); 
T     = I(1:s); 
x     = zeros(n,1); 
x(T)  = 0.1+abs(randn(s,1)); 
end