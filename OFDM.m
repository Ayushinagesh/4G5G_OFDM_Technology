close all;
clear all;

Nsub = 256;
Ncp = round(Nsub/4);
nBlocks = 10000;
L = 2;
EbdB = 1:4:55;
Eb = 10.^(EbdB/10);
No = 1;
SNR = 2 * Eb / No;
SNRdB = 10 * log10(SNR);
BER = zeros(size(EbdB));
BERt = zeros(size(EbdB));

for K = 1:length(SNRdB) 
    for blk = 1:nBlocks
        ChNoise = sqrt(No/2) * (randn(1, L + Nsub + Ncp - 1) + 1j * randn(1, L + Nsub + Ncp - 1));
        Bitsl = randi([0, 1], [1, Nsub]);
        BitsQ = randi([0, 1], [1, Nsub]);  
        Sym = (2 * Bitsl - 1) + 1j * (2 * BitsQ - 1);  
        h = 1/sqrt(2) * (randn(1, L) + 1j * randn(1, L));
        hFFT = fft(h, Nsub);

        LoadedSym = sqrt(Eb(K)) * Sym;  
        TxSamples = ifft(LoadedSym);
        TxSamCp = [TxSamples(Nsub-Ncp+1:Nsub),TxSamples];
        RxSamCP = conv(h,TxSamCp) + ChNoise;  

        RxSamples = RxSamCP(Ncp + 1 : Ncp + Nsub);
        R_Sym = fft(RxSamples, Nsub);
        Fout = R_Sym ./ hFFT;

        DecBits = real(Fout) > 0;
        DecBitsQ = imag(Fout) > 0;

        BER(K) = BER(K) + sum(DecBits ~= Bitsl) + sum(DecBitsQ ~= BitsQ);

    end
    BER(K) = BER(K) / (nBlocks * Nsub * 2);
end

SNReff = SNR * L / Nsub;  
BERt = 0.5 ./ SNReff;

semilogy(SNRdB, BER, 'g-', 'linewidth', 3.0, 'MarkerFaceColor', 'g', 'MarkerSize', 9.0);
hold on;
semilogy(SNRdB, BERt, 'ro', 'linewidth', 3.0, 'MarkerFaceColor', 'r', 'MarkerSize', 9.0);
grid on;
axis tight;
title('BER for OFDM');  
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
legend('Simulation', 'Theory');
