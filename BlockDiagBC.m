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

end

