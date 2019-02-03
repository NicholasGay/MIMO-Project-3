function [ Rsum ] = BlockDiagBCEqualPower( H,C,Ptx )
% function [ Rsum ] = BlockDiagBCEqualPower( H,C,Ptx )
%
% The function computes the maximum achievable rate of a MIMO BC with block
% diagonalization and equal power allocation
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% C: K x 1 cell array with noise covariance C_nk per cell
% Ptx: joint available transmit power
% Outputs
% Rsum: achievable sum rate

[M,N] = size(H{1});
K = length(H);
Rsumk = zeros(1,K);
for ik = 1:K %loop for users
    H_bar = [];
    %Front
    for ii = 1:ik-1
        H_bar = [H_bar, H{ii}'];
    end
    %Back
    for ii = ik+1:K
        H_bar = [H_bar,H{ii}'];
    end
    H_bar = H_bar';
    %% Finding Pk
    inverse_arg = H_bar*H_bar';
    Pk = eye(N) - H_bar'/inverse_arg*H_bar;
    
    %% Finding V
    [Vk,eigenvalue] = eig(Pk);
    eigenvalue = diag(eigenvalue);
    Vk(:,(abs(eigenvalue) <= 1e-6)) = [];
    %% H_hat
    H_hat = H{ik}*Vk;
    
    %% finding Tk
    [V,phi] = eig(H_hat'/C{ik}*H_hat);
    phi = diag(phi);
    [psi,~,~] = waterfilling(phi,Ptx/K);
    psi = diag(psi);
    %Tk = V*sqrtm(psi);
    Q = V*psi*V';
    %% finding rate
    second_argument = H_hat'/C{ik}*H_hat*Q;
    argument= eye(M) + second_argument;
    Rsumk(ik) = log(det(argument));
 end
Rsum = sum(Rsumk);
