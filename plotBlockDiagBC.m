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
RsumEqual = zeros(1,no_P);
RsumBlock = zeros(1,no_P);
for no = 1:no_P
  RsumEqual(no) = real(BlockDiagBCEqualPower( H,C,P(no) ));
  RsumBlock(no) = real(BlockDiagBC( H,C,P(no) ));
end

figure;
hold on;

plot(Ptx,RsumEqual,'DisplayName','RsumEqual');
plot(Ptx,RsumBlock ,'DisplayName','RsumBlock');
     
hold off;
xlabel('Ptx in [dB]');
ylabel('R in [bits/channel use]');
legend('show','Location','NorthWest');
