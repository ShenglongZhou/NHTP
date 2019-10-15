% This code presents the success recovery rate of NHTP
clc; clear; close all; 

test    = 1; %=1 succ rate v.s. s; =2 succ rate v.s. m/n
ExMat   = 1; %=1 Gaussian matrix;  =2 Partial DCT matrix

n       = 256; 
m       = ceil(0.25*n);
s       = ceil(0.05*n);
noS     = 500;
MatType = {'GaussianMat','PartialDCTMat'};
if test == 1; test0 = 6:2:36; 
else;         test0 = 0.08:0.02:0.22;
end
 
SuccRate     = [];
pars.display = 0;
pars.draw    = 0;
par.eta      = 2; 
for j    = 1:length(test0) 
    rate = 0; 
    if test==1; s = test0(j);
    else;       m = floor(test0(j)*n);
    end    
    for S = 1:noS         
        data = compressed_sensing_data(MatType{ExMat},m,n,s,0 );       
        func = @(x,fgh,T1,T2)compressed_sensing(x,fgh,T1,T2,data);
        out  = NHTP(n,s,func,pars); 
        clc; SuccRate     
        rate = rate+ (norm(out.sol-data.x_opt)/norm(data.x_opt)<1e-2); 
    end
    clc; SuccRate  = [SuccRate rate]    
end

figure
plot(test0,SuccRate/noS,'r*-') ; hold on
if test==1; xlabel('s')  ;
else;       xlabel('m/n');
end
ylabel('Success Rate');
axis([min(test0) max(test0) 0 1]); grid on;
legend('NHTP','Location','NorthEast'); hold on 
saveas(figure(1), 'outputs\SuccessRate.eps','epsc');
saveas(figure(1), 'outputs\SuccessRate.fig');
