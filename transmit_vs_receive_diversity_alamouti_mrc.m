% ============================================================
% Project Title:
% Transmit vs Receive Diversity (Alamouti vs MRC)
%
% Systems Compared:
% 1) No Diversity (1Tx,1Rx)
% 2) Transmit Diversity using Alamouti (2Tx,1Rx)
% 3) Receive Diversity using MRC (1Tx,2Rx)
%
% Modulation: BPSK
% Channel: Rayleigh Fading + AWGN
% Output: BER vs Eb/No
% ============================================================

clc;
clear;
close all;

frmLen = 100;       % frame length
numPackets = 1000;  % number of packets
EbNo = 0:2:20;      % Eb/No range

N = 2;              % number of Tx antennas
M = 2;              % number of Rx antennas

P = 2;              % BPSK modulation order

%% Modulator / Demodulator
hMod = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator('OutputDataType','double');

%% Alamouti Encoder / Combiner
hAlamoutiEnc = comm.OSTBCEncoder;
hAlamoutiDec = comm.OSTBCCombiner;

%% AWGN Channel
Hawgn1Rx = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
                            'SignalPower',1);

Hawgn2Rx = clone(Hawgn1Rx);

%% Error Calculators
hErrorCalc1 = comm.ErrorRate;
hErrorCalc2 = comm.ErrorRate;
hErrorCalc3 = comm.ErrorRate;

%% Fix random stream (corrected)
S = RandStream.create('mt19937ar','seed',55408);
prevStream = RandStream.setGlobalStream(S);

%% Channel Matrix
H = zeros(frmLen,N,M);

%% BER storage
ber_noDiver  = zeros(3,length(EbNo));
ber_Alamouti = zeros(3,length(EbNo));
ber_MaxRatio = zeros(3,length(EbNo));
ber_thy2     = zeros(1,length(EbNo));

%% Figure setup
figure
grid on
hold on
set(gca,'YScale','log')
xlabel('Eb/No (dB)')
ylabel('BER')
title('Transmit vs Receive Diversity')

%% Simulation loop
for idx = 1:length(EbNo)

    reset(hErrorCalc1)
    reset(hErrorCalc2)
    reset(hErrorCalc3)

    Hawgn1Rx.EbNo = EbNo(idx);
    Hawgn2Rx.EbNo = EbNo(idx);

    for packetIdx = 1:numPackets

        % Generate data
        data = randi([0 P-1], frmLen, 1);

        % BPSK modulation
        modData = step(hMod,data);

        % Alamouti encoding
        encData = step(hAlamoutiEnc,modData);

        % Rayleigh channel
        H(1:N:end,:,:) = (randn(frmLen/2,N,M) + ...
                          1i*randn(frmLen/2,N,M))/sqrt(2);

        H(2:N:end,:,:) = H(1:N:end,:,:);

        % Extract channels
        H11 = H(:,1,1);
        H21 = H(:,:,1)/sqrt(2);
        H12 = squeeze(H(:,1,:));

        % Channel outputs
        chanOut11 = H11 .* modData;
        chanOut21 = sum(H21 .* encData,2);
        chanOut12 = H12 .* repmat(modData,1,2);

        % Add noise
        rxSig11 = step(Hawgn1Rx,chanOut11);
        rxSig21 = step(Hawgn1Rx,chanOut21);
        rxSig12 = step(Hawgn2Rx,chanOut12);

        % Alamouti decoding
        decData = step(hAlamoutiDec,rxSig21,H21);

        % Detection
        demod11 = step(hDemod,rxSig11 .* conj(H11));
        demod21 = step(hDemod,decData);
        demod12 = step(hDemod,sum(rxSig12 .* conj(H12),2));

        % BER update
        ber_noDiver(:,idx)  = step(hErrorCalc1,data,demod11);
        ber_Alamouti(:,idx) = step(hErrorCalc2,data,demod21);
        ber_MaxRatio(:,idx) = step(hErrorCalc3,data,demod12);

    end

    % theoretical BER
    ber_thy2(idx) = berfading(EbNo(idx),'psk',2,2);

    % Plot curves
    semilogy(EbNo(1:idx),ber_noDiver(1,1:idx),'r-*','LineWidth',2)
    semilogy(EbNo(1:idx),ber_Alamouti(1,1:idx),'g-o','LineWidth',2)
    semilogy(EbNo(1:idx),ber_MaxRatio(1,1:idx),'b-s','LineWidth',2)
    semilogy(EbNo(1:idx),ber_thy2(1:idx),'m','LineWidth',2)

    legend('No Diversity (1Tx,1Rx)',...
           'Alamouti (2Tx,1Rx)',...
           'MRC (1Tx,2Rx)',...
           'Theoretical Diversity')

    drawnow

end

%% Restore random stream
RandStream.setGlobalStream(prevStream);