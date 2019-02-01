%% Question 8
clear all
close all
clc
load('exampleMIMOBCs.mat')

Ptx_dB = 10;
Ptx = 10^(Ptx_dB/10);

[Q,Csum] = DualMACSumRateMaximization(H,Ptx);
K = length(H);
order = 1:K;
S = MACtoBCtransform(Q,H,order);

BC_power = zeros(1,K);
for k = 1:K
    BC_power(k) = 10*log10(trace(S{k}));
end