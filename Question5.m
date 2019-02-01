clc
clear all
close all
%% Question 5
load('exampleMIMOBCs')
K = length(H);
order = K:-1:1;
S = MACtoBCtransform(Q,H,order);

[ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order );
