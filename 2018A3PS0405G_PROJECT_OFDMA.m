% Code to simulate SER vs SNR [dB] for OFDMA
clc
clear
SP.subband = 0;
SP.inputBlockSize = 32; % Input data block size
SP.FFTsize = 512; % The size of the FFT and IFFT.
SP.CPsize = 20; % Cyclic Prefic (CP) length.
SP.SNR = 0:2:30; % Simulated SNR range is from 0 dB to 30 dB.
SP.numRun = 10^4; % The number of simulation iterations is 10^4.
% Channels based on 3GPP TS 25.104.
pedAchannel = [1 10^(-9.7/20) 10^(-22.8/20)];
pedAchannel = pedAchannel/sqrt(sum(pedAchannel.^2)); % Normalize the channel.
vehAchannel = [1 0 10^(-1/20) 0 10^(-9/20) 10^(-10/20) 0 0 0 10^(-15/20) 0 0 0 10^(-20/20)];
vehAchannel = vehAchannel/sqrt(sum(vehAchannel.^2)); % Normalize the channel.
idenChannel = 1; % This is identity channel for AWGN simulation.
% Set the type of the channel.
SP.channel = vehAchannel;

% Run the simulation for OFDMA.
M = 64; % M-ary OAM
numSymbols = SP.FFTsize;
H_channel = fft(SP.channel,SP.FFTsize);
SER = zeros(1,length(SP.SNR));
totalSubcarriers=512;
Q = numSymbols/SP.inputBlockSize; % Bandwidth spreading factor of IFDMA.
Q_tilda = Q-5; % Bandwidth spreading factor of DFDMA. Q_tilda < Q.
for n=1:length(SP.SNR)
    errCount_lfdma = 0;
    errCount_ifdma = 0;
    errCount_dfdma = 0;
    for k=1:SP.numRun
        % Generate random data.
        inputData = randi([0 M-1], 1, SP.inputBlockSize);
        inputSymbols = qammod(inputData, M, 'UnitAveragePower', true);
        % Initialize the subcarriers.
        Y = zeros(1,numSymbols); %For LFDMA
        Y1 = zeros(1,numSymbols); %For IFDMA
        Y2 = zeros(1,numSymbols); %For DFDMA
        % Subcarrier mapping. 
        Y(1:SP.inputBlockSize) = inputSymbols; %LFDMA
        Y1(1:Q:numSymbols) = inputSymbols; %IFDMA
        Y2(1:Q_tilda:Q_tilda*SP.inputBlockSize) = inputSymbols; %DFDMA
        % Perform OFDM modulation using IFFT.
        TxSamples_lfdma = sqrt(SP.FFTsize)*ifft(Y);
        TxSamples_ifdma = sqrt(SP.FFTsize)*ifft(Y1);
        TxSamples_dfdma = sqrt(SP.FFTsize)*ifft(Y2);
        % Add CP.
        ofdmSymbol_lfdma = [TxSamples_lfdma(numSymbols-SP.CPsize+1:numSymbols) TxSamples_lfdma];
        ofdmSymbol_ifdma = [TxSamples_ifdma(numSymbols-SP.CPsize+1:numSymbols) TxSamples_ifdma];
        ofdmSymbol_dfdma = [TxSamples_dfdma(numSymbols-SP.CPsize+1:numSymbols) TxSamples_dfdma];
        % Propagate through multi-path channel.
        RxSamples_lfdma = filter(SP.channel, 1, ofdmSymbol_lfdma);
        RxSamples_ifdma = filter(SP.channel, 1, ofdmSymbol_ifdma);
        RxSamples_dfdma = filter(SP.channel, 1, ofdmSymbol_dfdma);
        % Generate AWGN with appropriate noise power.
        tmp = randn(2, numSymbols+1*SP.CPsize);
        complexNoise = (tmp(1,:) + 1i*tmp(2,:))/sqrt(2);
        noisePower = 10^(-SP.SNR(n)/10);
        % Add AWGN to the transmitted signal.
        RxSamples_lfdma = RxSamples_lfdma + sqrt(noisePower)*complexNoise;
        RxSamples_ifdma = RxSamples_ifdma + sqrt(noisePower)*complexNoise;
        RxSamples_dfdma = RxSamples_dfdma + sqrt(noisePower)*complexNoise;
        % Remove CP.
        EstSymbols_lfdma = RxSamples_lfdma(SP.CPsize+1:numSymbols+SP.CPsize);
        EstSymbols_ifdma = RxSamples_ifdma(SP.CPsize+1:numSymbols+SP.CPsize);
        EstSymbols_dfdma = RxSamples_dfdma(SP.CPsize+1:numSymbols+SP.CPsize);
        % Convert the received signal into frequency domain.
        Y_lfdma = fft(EstSymbols_lfdma, SP.FFTsize)/sqrt(SP.FFTsize);
        Y_ifdma = fft(EstSymbols_ifdma, SP.FFTsize)/sqrt(SP.FFTsize);
        Y_dfdma = fft(EstSymbols_dfdma, SP.FFTsize)/sqrt(SP.FFTsize);
        %Sub-carrier Demapping
        Y_lfdma = Y_lfdma(1:SP.inputBlockSize);
        Y_ifdma = Y_ifdma(1:Q:numSymbols);
        Y_dfdma = Y_dfdma(1:Q_tilda:Q_tilda*SP.inputBlockSize);
        % Perform channel equalization in the frequencydomain.
        H_eff = H_channel([1:SP.inputBlockSize]+SP.inputBlockSize*SP.subband);
        Y1_lfdma = Y_lfdma./H_eff;
        H_eff = H_channel(1+SP.subband:Q:numSymbols);
        Y1_ifdma = Y_ifdma./H_eff;
        H_eff = H_channel(1+SP.subband:Q_tilda:Q_tilda*SP.inputBlockSize);
        Y1_dfdma = Y_dfdma./H_eff;
        % Perform demodulation
        EstData_lfdma = qamdemod(Y1_lfdma, M, 'UnitAveragePower', true);
        EstData_ifdma = qamdemod(Y1_ifdma, M, 'UnitAveragePower', true);
        EstData_dfdma = qamdemod(Y1_dfdma, M, 'UnitAveragePower', true);
        % Error counting
        errCount_lfdma = errCount_lfdma + symerr(EstData_lfdma,inputData);
        errCount_ifdma = errCount_ifdma + symerr(EstData_ifdma,inputData);
        errCount_dfdma = errCount_ifdma + symerr(EstData_dfdma,inputData);
    end
    % Symbol error rate
    SER_lfdma(n) = errCount_lfdma / (numSymbols*SP.numRun);
    SER_ifdma(n) = errCount_ifdma / (numSymbols*SP.numRun);
    SER_dfdma(n) = errCount_ifdma / (numSymbols*SP.numRun);
end

% Plot results
figure
semilogy(SP.SNR,SER_lfdma,'r')
hold on
semilogy(SP.SNR,SER_ifdma,'b')
semilogy(SP.SNR,SER_dfdma,'g')
legend('LFDMA','IFDMA','DFDMA')
hold off
grid on
xlabel('SNR [dB]')
ylabel('SER')
title('OFDMA for AWGN Channel')

% Uncomment below code to get plot in the subplot format   
%   subplot(3,1,1)
%   semilogy(SP.SNR,SER_lfdma,'r')
%   xlabel('SNR [dB]')
%   ylabel('SER')
%   title('OFDM LFDMA')
%   
%   subplot(3,1,2)
%   semilogy(SP.SNR,SER_ifdma,'b')
%   xlabel('SNR [dB]')
%   ylabel('SER')
%   title('OFDM IFDMA')
%   
%   subplot(3,1,3)
%   semilogy(SP.SNR,SER_dfdma,'g')
%   xlabel('SNR [dB]')
%   ylabel('SER')
%   title('OFDM DFDMA')
