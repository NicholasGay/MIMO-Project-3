%% Question 7
clear all
close all
clc
load('exampleMIMOBCs.mat')

Ptx_dB = 10;
Ptx = 10^(Ptx_dB/10);

[Q,Csum] = DualMACSumRateMaximization(H,Ptx);