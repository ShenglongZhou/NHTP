function Out = NHTP(func,n,s,pars)

% This code aims at solving the sparsity constrained optimization with form
%
%         min_{x\in R^n} f(x)  s.t.  \|x\|_0<=s
%
% where s is the given sparsity, which is << n.  
%--------------------------------------------------------------------------
% Inputs:
%     func:  A function handle defines (objective,gradient,sub-Hessain) (required)
%     n   :  Dimension of the solution x                                (required)
%     s   :  Sparsity level of x, an integer between 1 and n-1          (required) 
%
%     pars:  Parameters are all OPTIONAL
%            pars.x0      --  Starting point of x,  pars.x0=zeros(n,1) (default)
%            pars.eta     --  A positive parameter, a default one is given related to inputs  
%            pars.display --  Display results or not for each iteration (default, 1)
%            pars.draw    --  Draw or not draw a graph (default, 0) 
%            pars.maxit   --  Maximum number of iterations, (default,2000) 
%            pars.tol     --  Tolerance of the halting condition, (default,1e-6)
%
% Outputs:
%     Out.sol :   The sparse solution x
%     Out.time:   CPU time
%     Out.iter:   Number of iterations
%     Out.obj :   Objective function value at Out.sol 
%--------------------------------------------------------------------------
% This code is programmed based on the algorithm proposed in 
% "S. Zhou, N. Xiu and H. Qi, Global and Quadratic Convergence of Newton 
% Hard-Thresholding Pursuit, Journal of Machine Learning Research, 2019."
% Send your comments and suggestions to <<< slzhou2021@163.com >>> 
% Warning: Accuracy may not be guaranteed !!!!! 
%--------------------------------------------------------------------------

warning off;
t0  = tic;
if  nargin<3
    disp(' No enough inputs. No problems will be solverd!'); return;
end
if nargin < 4; pars = [];  end 
if isfield(pars,'display');display = pars.display;else; display = 1;    end
if isfield(pars,'draw');   draw    = pars.draw;   else; draw    = 0;    end
if isfield(pars,'maxit');  maxit   = pars.maxit;  else; maxit   = 2000; end
if isfield(pars,'tol');    tol     = pars.tol;    else; tol     = 1e-6; end  
if isfield(pars,'x0');     x0      = pars.x0;     else; x0 = zeros(n,1);end 
 
[x0,obj,g,eta] = getparameters(n,s,x0,func,pars);  
x       = x0;
beta    = 0.5;
sigma   = 5e-5;
delta   = 1e-10;
pcgtol  = 0.1*tol*s;
T0      = [];
Error   = zeros(1,maxit);
Obj     = zeros(1,maxit);
FNorm   = @(var)norm(var,'fro')^2;
xo      = zeros(n,1);

if  display 
    fprintf(' Start to run the solver -- NHTP \n');
    fprintf(' ------------------------------------------------\n');
    fprintf(' Iter       Error         Objective        Time \n'); 
    fprintf(' ------------------------------------------------\n');
end

% Initial check for the starting point
if  FNorm(g)<1e-20 && nnz(x)<=s
    fprintf(' Starting point is a good solution. Stop NHTP\n'); 
    Out.sol  = x;
    Out.obj  = obj;
    Out.time = toc(t0);
    return;
end

if  max(isnan(g))
    x0      = zeros(n,1);
    rind    = randi(n);
    x0(rind)= rand;
    [obj,g] = func(x0,[],[]);
end
 
% The main body  
for iter = 1:maxit
     
    xtg   = x0-eta*g;
    [~,T] = maxk(abs(xtg),s);  
    T     = sort(T);
    TTc   = setdiff(T0,T);
    flag  = isempty(TTc);    
    gT    = g(T);
    
    % Calculate the error for stopping criteria   
    xtaus           = max(0,max(abs(g))-min(abs(x(T)))/eta);
    if  flag
        FxT         = sqrt(FNorm(gT));
        Error(iter) = xtaus + FxT;
    else
        FxT         = sqrt(FNorm(gT)+ abs(FNorm(x)-FNorm(x(T))) );
        Error(iter) = xtaus + FxT;    
    end
     
    if  display 
        fprintf('%4d       %5.2e       %5.2e      %6.3fsec\n',...
                iter,Error(iter),obj,toc(t0)); 
    end
             
    % Stopping criteria
    if Error(iter)<tol; break;  end
    
    % update next iterate
    if  iter   == 1 || flag           % update next iterate if T==supp(x^k)   c
        H       =  func(x0,T,[]); 
        if ~isa(H,'function_handle')
            d   = H\(-gT);
        else
            d   = my_cg(H,-gT,pcgtol,50,zeros(s,1)); 
        end
        dg      = sum(d.*gT);
        ngT     = FNorm(gT);
        if dg   > max(-delta*FNorm(d), -ngT) || isnan(dg) 
            d   = -gT; 
            dg  = ngT; 
        end
    else                              % update next iterate if T~=supp(x^k) 
        [H,D]   = func(x0,T,TTc);
        if isa(D,'function_handle')
           Dx  = D(x0(TTc)); 
        else       
           Dx  = D*x0(TTc); 
        end        
        if isa(H,'function_handle')
           d    = my_cg(H,Dx-gT,pcgtol,50,zeros(s,1)); 
        else       
           d    = H\(Dx-gT); 
        end              
        Fnz     = FNorm(x(TTc))/4/eta;
        dgT     = sum(d.*gT);
        dg      = dgT-sum(x0(TTc).*g(TTc));
        
        delta0  = delta;
        if Fnz  > 1e-4; delta0 = 1e-4; end
        ngT     = FNorm(gT);
        if dgT  > max(-delta0*FNorm(d)+Fnz, -ngT) || isnan(dg) 
            d   = -gT; 
            dg  = ngT; 
        end
    end
    
    alpha    = 1; 
    x        = xo;    
    obj0     = obj;        
    Obj(iter)= obj;
    
    % Amijio line search

    for i      = 1:6
        x(T)   = x0(T) + alpha*d;
        obj    = func(x,[],[]);
        if obj < obj0  + alpha*sigma*dg; break; end        
        alpha  = beta*alpha;
    end
    
    % Hard Thresholding Pursuit if the obj increases
    fhtp    = 0;
    if obj  > obj0 
       x(T) = xtg(T); 
       obj  = func(x,[],[]); 
       fhtp = 1;
    end
    
    % Stopping criteria
    flag1   = (abs(obj-obj0)<1e-6*(1+abs(obj)) && fhtp); 
    flag2   = (abs(obj-obj0)<1e-10*(1+abs(obj))&& Error(iter)<1e-2);
    if  iter>10 &&  (flag1 || flag2)      
        if obj > obj0
           iter    = iter-1; 
           x       = x0; 
           T       = T0; 
        end   
        break;
     end 
 
    T0      = T; 
    x0      = x; 
    [obj,g] = func(x,[],[]);
    
    % Update eta
    if  mod(iter,50)==0  
        if Error(iter)>1/iter^2  
        if iter<1500; eta = eta/1.05; 
        else;         eta = eta/1.5; 
        end     
        else;         eta = eta*1.25;   
        end
    end     
end


% results output
Out.time    = toc(t0);
Out.iter    = iter;
Out.sol     = x;
Out.obj     = obj; 

if  draw
    figure
    subplot(121) 
    Obj(iter)= obj;  
    PlotFun(Obj,iter,'r.-','f(x^k)'); 
    subplot(122) 
    PlotFun(Error,iter,'r.-','error') 
end


if display 
   normgrad = sqrt(FNorm(g)); 
   fprintf(' ------------------------------------------------\n');
   if normgrad<1e-5
      fprintf(' A global optimal solution might be found\n');
      fprintf(' because of ||gradient|| = %5.2e!\n', normgrad); 
      if Out.iter>1500
      fprintf('\n Since the number of iterations reaches to %d\n',Out.iter);
      fprintf(' Try to rerun the solver with setting a smaller pars.eta \n'); 
      end
      fprintf(' ------------------------------------------------\n');
   end
end


end

% initialize parameters ---------------------------------------------------
function [x0,obj,g,eta]=getparameters(n,s,x0,func,pars)

    if isfield(pars,'x0') && norm(x0)>0
       [obj0,g0] = func(zeros(n,1),[],[]);  
       [obj,g]   = func(pars.x0,[],[]); 
       if obj0   < obj/10
          x0     =  zeros(n,1); 
          obj    = obj0;  
          g      = g0; 
       else  
          ns0    = nnz(pars.x0);
          if  ns0==s
              [~,T]    = maxk(pars.x0,s,'ComparisonMethod','abs'); 
              x0       = pars.x0;  
              pars.eta = min(abs(x0(T)))/(1+max(abs(g(setdiff(1:n, T)))));   
          elseif ns0<s
              x0       = pars.x0;  
              pars.eta = max(x0(x0>0.1))/(1+max(abs(g)));   
          else 
              [~,T]    = maxk(pars.x0,s,'ComparisonMethod','abs'); 
              x0       = zeros(n,1);
              x0(T)    = pars.x0(T);  
              pars.eta = max(x0(x0>0.1))/(1+max(abs(g)));  
          end
          
          if isempty(pars.eta) 
          pars.eta = max(abs(x0))/(1+max(abs(g))); 
          end          
       end
    else
        [obj,g]  = func(x0,[],[]); 
    end
 
    
    if isfield(pars,'eta')      
        eta    = pars.eta;       
    else % set a proper parameter eta
        [~,g1] = func(ones(n,1),[],[]) ;
        abg1   = abs(g1);
        T      = find(abg1>1e-8);
        maxe   = sum(1./(abg1(T)+eps))/nnz(T);
        if  isempty(T) 
            eta = 10*(1+s/n)/min(10, log(n));
        else
            if maxe>2
                eta  = (log2(1+maxe)/log2(maxe))*exp((s/n)^(1/3));
            elseif maxe<1
                eta  = (log2(1+ maxe))*(n/s)^(1/2);    
            else
                eta  = (log2(1+ maxe))*exp((s/n)^(1/3));
            end     
        end
    end 
    
end

% plot the graph: iter v.s. obj (or error) --------------------------------
function  PlotFun(input,iter,c, key) 
    if  input(iter)>1e-40 && input(iter)<1e-5
        semilogy(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
    else
        plot(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
    end
    xlabel('Iter'); ylabel(key); grid on    
    
end

% conjugate gradient-------------------------------------------------------
function x = my_cg(fx,b,cgtol,cgit,x)
    if ~isa(fx,'function_handle'); fx = @(v)fx*v; end
    r = b;
    if nnz(x)>0; r = b - fx(x);  end
    e = norm(r,'fro')^2;
    t = e;
    p = r;
    for i = 1:cgit  
        if e < cgtol*t; break; end
        w  = fx(p);
        pw = p.*w;
        a  = e/sum(pw(:));
        x  = x + a * p;
        r  = r - a * w;
        e0 = e;
        e  = norm(r,'fro')^2;
        p  = r + (e/e0)*p;
    end 
end
