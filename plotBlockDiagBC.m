% The Script calculates the sum capacity of the MIMO BC and the maximum
% achievable rates with block diagonalization with equal power allocation
% and block diagonlization with optimzal power allocation and plots the
% result versus Ptx

load('exampleBlockDiagMIMOBC.mat')

% transmit powers
Ptx = -10:2:40;
P = 10.^(Ptx/10);
no_P = length(P);

% Rate calculations

for no = 1:no_P
  
end

figure;
hold on;


hold off;
xlabel('Ptx in [dB]');
ylabel('R in [bits/channel use]');
legend('show','Location','NorthWest');
