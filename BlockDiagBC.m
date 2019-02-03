function [ Rsum ] = BlockDiagBC( H,C,Ptx )
% function [ Rsum ] = BlockDiagBC( H,C,Ptx )
%
% The function computes the maximum achievable rate of a MIMO BC with block
% diagonalization.
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% C: K x 1 cell array with noise covariance C_nk per cell
% Ptx: joint available transmit power
% Outputs
% Rsum: achievable sum rate

[M,N] = size(H{1});
K = length(H);

%% FINDING V
V_vec = [];
H_vec = [];
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
    
    V_vec = [V_vec,Vk];
    
    H_vec = [H_vec;H{ik}];
    
end

C_block = blkdiag(C{:});

H_hat = H_vec*V_vec;

%% finding Tk
[V_t,phi] = eig(H_hat'/C_block*H_hat);
phi = diag(phi);
[psi,~,~] = waterfilling(phi,Ptx);
psi = diag(psi);
%T = V_t*sqrtm(psi);
Q = V_t*psi*V_t';
%% finding rate
second_argument = H_hat'/C_block*H_hat*Q;
argument= eye(M*K) + second_argument;
Rsum = log(det(argument));

end

