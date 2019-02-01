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


end
