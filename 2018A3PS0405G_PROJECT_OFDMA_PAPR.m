% Code to simuulate PAPR for OFDMA
clc
clear

M = 16; % M-ary QAM
totalSubcarriers = 512; % Number of total subcarriers.
numSymbols = 32; % Input Data block size.
Q = totalSubcarriers/numSymbols; % Bandwidth spreading factor of IFDMA.
Q_tilda = 15; % Bandwidth spreading factor of DFDMA. Q_tilda < Q.
numRuns = 1e5; % Number of runs.
papr = zeros(1,numRuns); % Initialize the PAPR results. (LFDMA)
papr1 = zeros(1,numRuns); % Initialize the PAPR results. (IFDMA)
papr2 = zeros(1,numRuns); % Initialize the PAPR results. (DFDMA)
for n = 1:numRuns
    % Generate random data.
    inputData = randi([0 M-1], 1, numSymbols);
    inputSymbols = qammod(inputData, M, 'UnitAveragePower', false);
    % Initialize the subcarriers.
    Y = zeros(totalSubcarriers,1); %For LFDMA
    Y1 = zeros(totalSubcarriers,1); %For IFDMA
    Y2 = zeros(totalSubcarriers,1); %For DFDMA
    % Subcarrier mapping. LFDMA
    Y(1:numSymbols) = inputSymbols;
    % Subcarrier mapping. IFDMA
    Y1(1:Q:totalSubcarriers) = inputSymbols;
    % Subcarrier mapping. DFDMA
    Y2(1:Q_tilda:Q_tilda*numSymbols) = inputSymbols;
    % IDFT
    y = ifft(Y);
    y1 = ifft(Y1);
    y2 = ifft(Y2);
    % Calculate PAPR.
    papr(n) = 10*log10(max(abs(y).^2)/mean(abs(y).^2));
    papr1(n) = 10*log10(max(abs(y1).^2)/mean(abs(y1).^2));
    papr2(n) = 10*log10(max(abs(y2).^2)/mean(abs(y2).^2));
end

% Plot CCDF.
[N,X] = histcounts(papr, 200,'Normalization','count');
[N1,X1] = histcounts(papr1, 200,'Normalization','count');
[N2,X2] = histcounts(papr2, 200,'Normalization','count');
figure
 semilogy(X(1:end-1),1-cumsum(N)/max(cumsum(N)),'-r') % For LFDMA
 hold on
 semilogy(X1(1:end-1),1-cumsum(N1)/max(cumsum(N1)),'-g') % For IFDMA
 hold on
 semilogy(X2(1:end-1),1-cumsum(N2)/max(cumsum(N2)),'-b') % For DFDMA
 hold off
 legend('LFDMA','IFDMA','DFDMA')
 grid on
 xlabel('PAPR')
 ylabel('CCDF')
 title('OFDMA PAPR')

% Uncomment below code to get plot in the subplot format 
%  subplot(3,1,1)
%  semilogy(X(1:end-1),1-cumsum(N)/max(cumsum(N)),'-r') % For LFDMA
%  xlabel('PAPR')
%  ylabel('CCDF')
%  title('LFDMA - OFDMA PAPR')
%  subplot(3,1,2)
%  semilogy(X1(1:end-1),1-cumsum(N1)/max(cumsum(N1)),'-g') % For IFDMA
%  xlabel('PAPR')
%  ylabel('CCDF')
%  title('IFDMA - OFDMA PAPR')
%  subplot(3,1,3)
%  semilogy(X2(1:end-1),1-cumsum(N2)/max(cumsum(N2)),'-b') % For DFDMA
%  xlabel('PAPR')
%  ylabel('CCDF')
%  title('DFDMA - OFDMA PAPR')