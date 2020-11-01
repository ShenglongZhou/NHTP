% This code demonstrates the success recovery rate in compressed sensing
clc; clear; close all; 

test    = 1; %=1 succ rate v.s. s; =2 succ rate v.s. m/n
ExMat   = 1; %=1 Gaussian matrix;  =2 Partial DCT matrix

n       = 256; 
m       = ceil(0.25*n);
s       = ceil(0.05*n);
noS     = 100;
MatType = {'GaussianMat','PartialDCTMat'};
switch  test
case 1; sm = ceil(linspace(6,36,15));
case 2; sm = linspace(0.1,0.3,12);
end    
    
SucRate      = [];
pars.display = 0;
pars.draw    = 0;
par.eta      = 5; 
for j        = 1:length(sm)
    rate     = 0; 
    switch  test
    case 1; s = sm(j);
    case 2; m = ceil(sm(j)*n);
    end    
    for S = 1:noS         
        data = compressed_sensing_data(MatType{ExMat},m,n,s,0 );       
        func = @(x,fgh,T1,T2)compressed_sensing(x,fgh,T1,T2,data);
        out  = NHTP(n,s,func,pars); clc; SucRate     
        rate = rate + (norm(out.sol-data.x_opt)/norm(data.x_opt)<1e-2); 
    end
    clc; SucRate  = [SucRate rate]    
end

xlab = {'s','m/n'};
figure(1), plot(sm,SucRate/noS,'r*-'), 
xlabel(xlab{test}), ylabel('Success Rate') 
axis([min(sm) max(sm) 0 1]); grid on;
legend('NHTP','Location','NorthEast'); hold on 
saveas(figure(1), 'outputs\SuccessRate.eps','epsc');
saveas(figure(1), 'outputs\SuccessRate.fig');
