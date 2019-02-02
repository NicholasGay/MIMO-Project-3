function [ Rsato ] = sato( H, Ptx )
% function [ Rsato ] = sato( H, Ptx )
%
% The function computes the Sato bound for a given channel H and transmit 
% power Ptx.
%
% Inputs
% H: M x N x 2 array of given channels to the users
% Ptx: Transmit power
%
% Outputs
% Rsato: Sato Bound

% Set YALMIP settings
options = sdpsettings('solver','sdpt3','verbose',0);

K = length(H);
[M,N] = size(H{1});

H_hat = [];
 for i = 1:K
     H_hat = vertcat(H_hat,H{i});
 end

 %% Initialize optimization variables
Z = sdpvar(M*K,M*K,'hermitian','complex');
K_var = sdpvar(N,N,'hermitian','complex');
miu = sdpvar(1,1,'hermitian','complex');

%% Define Constraint Set
Contraints = [Z>=0,K_var>=0,miu>=0,Z>=H_hat*K_var*H_hat'];

for i = 1:K
    S = [zeros(M,(i-1)*M),eye(M),zeros(M,(K-i)*M)];
    Contraints = [Contraints, S*Z*S.' <= miu*eye(M)];
end

%% Create Objective
Objective = (-logdet(K_var) + trace(K_var) + miu*Ptx - N);

%% Solve
sol = optimize(Contraints,Objective,options);

%% Retreve solution
K_var = value(K_var);
miu = value(miu);
Z = value(Z);

Rsato = (1/log(2))*(-log(det(K_var))+trace(K_var)+miu*Ptx-N);