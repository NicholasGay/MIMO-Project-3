function [psi,mu,K] = waterfilling(phi,Etx)
% Function [psi,mu] = waterfilling(phi,Etx)
%
% Waterfilling power allocation procedure under the total transmit
% power constraint sum(psi)<=Etx.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Etx: available sum transmit power Etx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% mu: value of the optimal waterlevel mu^prime
% K: number of active (non-zero) data streams K

phi = real(phi(:));
phi(phi<0) = 0;

% determine N
N = length(phi);

% sort phi in descending order
phi_sort = sort(phi,'descend');

% compute the array mu_arr for all K=1:N
mu_arr = (Etx + cumsum(1./phi_sort))./(1:N)';

% determine K
K = find(mu_arr > 1./phi_sort,1,'last');

% compute mu and psi
mu = mu_arr(K);
psi = max(zeros(N,1), mu - 1./phi);
