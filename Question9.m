%% Question 9 
clc
close all
clear all
%load example
load('exampleMIMOBCs')
Ptx = -15:5:30;
fig = plotSumRateBC(H,Ptx);

%% SLOPE
% Getting the Y data 
dataObjs = findobj(fig,'-property','YData');
y = dataObjs(1).YData;

slope = (y(10)-y(9))/5;
