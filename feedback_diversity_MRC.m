% Parameters
numBits = 1e7; % Number of bits to transmit
SNRdB = 0:2:20; % SNR range in dB
numDiversityBranches = 2; % Number of diversity branches

% Generate random bits
dataBits = randi([0 1], numBits, 1);

% BPSK Modulation
txSymbols = 2*dataBits - 1; % BPSK symbols: 0 -> -1, 1 -> 1

% Rayleigh Fading Channel
h = (randn(numBits, numDiversityBranches) + 1i*randn(numBits, numDiversityBranches)) / sqrt(2);

% Add AWGN Noise
noise = (randn(numBits, numDiversityBranches) + 1i*randn(numBits, numDiversityBranches)) / sqrt(2);

% Initialize BER array
BER = zeros(length(SNRdB), 1);

for idx = 1:length(SNRdB)

    SNR = 10^(SNRdB(idx)/10); % Convert SNR from dB to linear scale

    % Channel and noise scaling based on SNR
    rxSymbols = h .* repmat(txSymbols, 1, numDiversityBranches) + (1/sqrt(SNR)) * noise;

    % Calculate the combined SNR for MRC
    combinedSNR = sum(abs(h).^2, 2) * SNR;

    % Feedback power adjustment
    adjustedTxSymbols = txSymbols .* sqrt(combinedSNR);

    % Recompute received symbols
    rxSymbols = h .* repmat(adjustedTxSymbols, 1, numDiversityBranches) + (1/sqrt(SNR)) * noise;

    % Maximal Ratio Combining
    combinedSymbols = sum(conj(h) .* rxSymbols, 2) ./ sum(abs(h).^2, 2);

    % Decision
    rxBits = real(combinedSymbols) > 0;

    % BER calculation
    BER(idx) = sum(rxBits ~= dataBits) / numBits;

end

% Plot BER vs SNR
figure
semilogy(SNRdB, BER, 'b-o','LineWidth',2)
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate (BER)')
title('BER vs SNR for Feedback Diversity System')
legend('Maximal Ratio Combining with Feedback')