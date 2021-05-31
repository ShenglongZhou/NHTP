% This code demonstrates the success recovery rate in compressed sensing
clc; clear; close all; 

test    = 1; %=1 succ rate v.s. s; =2 succ rate v.s. m/n
ExMat   = 1; %=1 Gaussian matrix;  =2 Partial DCT matrix

n       = 256; 
m       = ceil(0.25*n);
s       = ceil(0.05*n);
noS     = 200;
MatType = {'GaussianMat','PartialDCTMat'};
switch  test
case 1; sm = ceil(linspace(6,36,16));
case 2; sm = linspace(0.1,0.3,12);
end    
    
SucRate      = [];
pars.display = 0;
pars.draw    = 0;
for j        = 1:length(sm)
    rate     = 0; 
    switch  test
    case 1; s = sm(j);
    case 2; m = ceil(sm(j)*n);
    end    
    for S = 1:noS         
        data = compressed_sensing_data(MatType{ExMat},m,n,s,0 );       
        func = @(x,T1,T2)CS(x,T1,T2,data);
        out  = NHTP(func,n,s,pars); clc; SucRate     
        rate = rate + (norm(out.sol-data.xopt)/norm(data.xopt)<1e-2); 
    end
    clc; SucRate  = [SucRate rate]  
    
    figure(1)
    set(gcf, 'Position', [1000, 200, 400 350]);
    xlab = {'s','m/n'};
    plot(sm(1:j),SucRate/noS,'r*-','LineWidth',1), 
    xlabel(xlab{test}), ylabel('Success Rate') 
    axis([min(sm) max(sm) 0 1]); grid on; 
    legend('NHTP','Location','NorthEast'); hold on, pause(0.1)
    
end


