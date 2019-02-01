function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
% function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
%
% The function calculates the sum capacity and the corresponding
% transmit covariance matrices of the K user MIMO MAC with a joint sum
% transmit power constraint.
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% Ptx: joint transmit power
%
% Outputs
% Q: K x 1 cell array of users' optimal transmit covariance Qk per cell
% Csum: sum capacity of the MIMO MAC
[M,N] = size(H{1});
K = length(H);

% Set YALMIP options
options = sdpsettings('solver','sdpt3','verbose',0);

% Initialize optimization variables
Q = cell(K,1);
for k = 1:K
    Q{k} = sdpvar(M,M,'hermitian','complex');
end

% Define Constraint Set
power = 0;
Constraints = zeros(M,M);
for k = 1:K
    power = power+trace(Q{k});
    Constraints = [Constraints, Q{k} >= 0];
end

Constraints = [Constraints,power<=Ptx];

% Objective
Z=eye(N);
for k = 1:K
    Z = Z+H{k}'*Q{k}*H{k};
end
Objective = -logdet(Z);

% Solve
optimize(Constraints, Objective,options);

% Solution
for k = 1:K
    Q{k} = value(Q{k});
end

%% Sum Capacity
sum_argument = zeros(N,N);
for k = 1:K
    sum_argument = sum_argument + H{k}'*Q{k}*H{k};
end
argument = eye(N)+sum_argument;
Csum = log2(det(argument));

end

