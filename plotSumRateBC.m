function [fig] = plotSumRateBC(H,Ptx,fig)
% Function [fig] = plotSumRateBC(H,Ptx,fig)
%
% The function plots the sum capacity Csum and the individual
% achievable rates of each user in the MIMO BC for the encoding orders
% (1,...,K) and (K,...,1).
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% Ptx: column vector of Ptx values in dB
% fig: figure handle (optional)
% Outputs
% fig: figure handle

% system dimensions
[M,N] = size(H{1});
K = length(H);

% transmit powers
Ptx = Ptx(:);
P = 10.^(Ptx/10);
no_P = length(P);

%% Programming Task 5------------------------------------------------------
%Init

S = cell(2,1);
order{1}=1:1:K;
order{2}=K:-1:1;
CsumBC = zeros(1,no_P);
CsumReversed = zeros(1,no_P);
Rate = zeros(3,no_P);
Rate_reversed = zeros(3,no_P);
%--------------------------------------------------------------------------


% Rate calculations
Csum = zeros(no_P,1);
for no = 1:no_P
    [Q,Csum(no)] = DualMACSumRateMaximization(H,P(no));
    %% Programming Task 5
    %----------------------------------------------------------------------
    for k = 1:2
        S{k} = MACtoBCtransform(Q,H,order{k});
        [RsumBC{k}, ~ ] = MAC_BC_rates( H, Q, S{k}, order{k} );
    end
    CsumBC(no) = sum(RsumBC{1});
    CsumReversed(no) = sum(RsumBC{2});
    Rate(:,no) = RsumBC{k};
    Rate_reversed(:,no) = RsumBC{k};
    %----------------------------------------------------------------------
    
end



if nargin<3, fig = figure; end
figure(fig);
hold on;
plot(Ptx,Csum, ...
     'b-','LineWidth',1.5, ...
     'DisplayName',['Csum']);
 
%% Pregramming Task 5------------------------------------------------------
plot(Ptx,CsumBC, ...
     'c-','LineWidth',1.5, ...
     'DisplayName',['Csum_BC']);

plot(Ptx,CsumReversed, ...
     'm-','LineWidth',1.5, ...
     'DisplayName',['Csum_BC_Reversed']);
 
%Individual Rates
%normal
plot(Ptx,Rate(1,:), ...
     'g:','LineWidth',1.5, ...
     'DisplayName',['User 1']);

plot(Ptx,Rate(2,:), ...
     'g--','LineWidth',1.5, ...
     'DisplayName',['User 2']);
 
plot(Ptx,Rate(3,:), ...
     'g-.','LineWidth',1.5, ...
     'DisplayName',['User 3']);
%Reversed
plot(Ptx,Rate_reversed(1,:), ...
     'k:','LineWidth',1.5, ...
     'DisplayName',['User 1_R']);

plot(Ptx,Rate_reversed(2,:), ...
     'k--','LineWidth',1.5, ...
     'DisplayName',['User 2_R']);
 
plot(Ptx,Rate_reversed(3,:), ...
     'k-.','LineWidth',1.5, ...
     'DisplayName',['User 3_R']);
%--------------------------------------------------------------------------
hold off;
xlabel('Ptx in [dB]');
ylabel('R in [bits/channel use]');
legend('show','Location','NorthWest');
