function ReoveryShow(xo,x,pos,ind)

    figure('Renderer', 'painters', 'Position', pos)
    axes('Position', [0.05 0.1 0.9 0.8] );
    stem(find(xo),xo(xo~=0),'bo-','MarkerSize',6, 'LineWidth',1),hold on
    stem(find(x),x(x~=0),'r*:', 'MarkerSize',4, 'LineWidth',1),hold on
    grid on, ymin= -0.1; ymax=0.2;
    xx  = [xo; x];
    if nnz(xx<0)>0, ymin= min(xx(xx<0))-0.1; end
    if nnz(xx>0)>0, ymax= max(xx(xx>0))+0.1; end   
    axis([1 length(x)  ymin  ymax])
    if ind
       snr   = norm(x-xo)/norm(x);
       st1   = strcat('Recovery accuracy = ',num2str(snr,4));         
       title(strcat(st1))
       set(0,'DefaultAxesTitleFontWeight','normal');

       legend('Ground-Truth', 'Recovered', 'Location', 'best')
    end
end

