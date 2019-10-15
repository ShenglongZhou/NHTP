function  PlotFun(input,iter,c, key) 
    if  input(iter)>1e-40 && input(iter)<1e-5
        semilogy(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
    else
        plot(1:iter,input(1:iter),c,'MarkerSize',7,'LineWidth',1);
    end
    xlabel('Iter'); ylabel(key); grid on    
    
end

