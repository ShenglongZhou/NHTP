% demon compressed sensing problems with real data
clc; clear; close all;

load 'DrivFace.mat';
load 'nlab.mat';   %'identity.mat';

[m,n]     = size(A);
s         = 10;
data.A    = A/sqrt(m); clear A
data.b    = y/sqrt(m); clear y
data.At   = data.A';
pars.eta  = 1;
func      = @(x,fgh,T1,T2)compressed_sensing(x,fgh,T1,T2,data);
out       = NHTP(n,s,func,pars)  

fprintf('\nSample size:       m=%4d,n=%4d\n', m,n);
fprintf('CPU time:         %6.3fsec\n',  out.time);
fprintf('Objective value:   %5.3e\n\n', out.obj);
saveas(figure(1), 'outputs\DriverFaces.fig');
saveas(figure(1), 'outputs\DriverFaces.eps','epsc');

