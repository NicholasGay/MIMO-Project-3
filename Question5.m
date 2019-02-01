clc
clear all
close all
%% Question 5
load('exampleMIMOBCs')
K = length(H);
%order = 1:1:K;
order = K:-1:1;
S = MACtoBCtransform(Q,H,order);

[ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order );
BC_power = zeros(1,K);
MAC_power = zeros(1,K);
for i = 1:K
    BC_power(i) = trace(real(S{i}));
    MAC_power(i) = trace(real(Q{i}));
end