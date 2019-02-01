function [S] = MACtoBCtransform(Q,H,order)
% Function [S] = MACtoBCtransform(Q,H,order)
% 
% Rate based duality transformation from the
% multiple access channel to the broadcast channel.
%
% Inputs
% Q: K x 1 cell array with matrix Qk per cell
% H: M x N cell array with channel Hk per cell
% order: 1 x K vector with BC encoding order
% Output
% S: K x 1 cell array of covariance matrix Sk per cell

% parameters
[M,N] = size(H{1});
K = length(H);

% initialization
S = cell(K,1);

%% Getting BC_order
[~,BC_order] = sort(order);

%% Starting loop
for i = K:-1:1
    user_k = BC_order(i);
    %FK
    sum_term = zeros(N,N);
    for ik = i+1:K
        sum_term = sum_term+S{BC_order(ik)};
    end
    Fk = H{user_k}*sum_term*H{user_k}'+eye(M);
    
    %Xk
    sum_term = zeros(N,N);
    for ik = 1:i-1
        sum_term = sum_term+H{BC_order(ik)}'*Q{BC_order(ik)}*H{BC_order(ik)};
    end
    Xk = eye(N)+sum_term;
    
    %Heff
    Heff =sqrtm(Fk)'\H{user_k}/sqrtm(Xk)';
    
    %Qkeff
    Qkeff = sqrtm(Fk)*Q{user_k}*sqrtm(Fk)';

    %Skeff
    Skeff = ptpTransform(Qkeff,Heff);
    
    S{user_k} = sqrtm(Xk)'\Skeff/sqrtm(Xk);
    
end

