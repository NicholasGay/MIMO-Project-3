function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
% function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
%
% The function calculates the achievable rates of the MAC and
% BC for given channels Hk and transmit covariance matrices Qk and Sk and
% a given decoding and encoding order.
%
% Input Specifications:
% H: Kx1 cell array with channel Hk per cell
% Q: Kx1 cell array with MxM matrix Qk per cell
% S: Kx1 cell array with NxN matrix Sk per cell
% order: 1xK vector in {1,...,K} with the BC encoding order
%
% Output Specifications:
% R_BC:  Kx1 array with achievable rates of the BC
% R_MAC: Kx1 array with achievable rates of the MAC

[M,N] = size(H{1});
K = length(H);

R_BC = zeros(K,1);
R_MAC = zeros(K,1);

%% Sorting
[~, indices_BC] = sort(order,'ascend');
indices_MAC = fliplr(indices_BC);

%% BC ---------------------------------------------------------------------
for i = 1:K
    k = indices_BC(i); %so i dont have to keep writing indices_BC
    
    %Sum Term
    sum_term = zeros(N,N);
    for iSum_term = i+1:K
        sum_term = sum_term+S{indices_BC(iSum_term)};
    end
    
    %interference term
    interference_term = H{k}*sum_term*H{k}'+eye(M);
    second_argument = H{k}'/(interference_term)*H{k}*S{k};
    R_BC(k,1) = log2(det(eye(N)+second_argument));
end



%% ENDOF BC ---------------------------------------------------------------

%% MAC---------------------------------------------------------------------


for ii = 1:K
    k = indices_MAC(ii);
    
    %Sum term
    sum_term = zeros(N,N);
    for iSum_term = ii+1:K
        i = indices_MAC(iSum_term);
        sum_term = sum_term+H{i}'*Q{i}*H{i};
    end
    
    %interference term
    interference_term = sum_term+eye(N);
    second_argument = H{k}/(interference_term)*H{k}'*Q{k};
    R_MAC(k,1) = log2(det(eye(M)+second_argument));

end







