clc;
clear;
close all;

%% =========================
%  MRC CODE
% =========================

N=10^3;
data=randi([0,1],1,N);
x=2*data-1;

nRx_max=20;
nRx=1:nRx_max;

snr_dB=1:5;

snr_sim_MRC=zeros([length(snr_dB) nRx_max]);

for j=1:nRx_max
    for k=1:length(snr_dB)

        h=randn(j,N)+(randn(j,N)*1i);
        x_kron=kron(ones(nRx(j),1),x);
        c_in=h.*x_kron;

        y=awgn(c_in,snr_dB(k),'measured');

        y_rec=sum(conj(h).*y,1);

        snr_sim_MRC(k,j)=mean(abs(y_rec));

    end
end

snr_sim_log_MRC=10*log10(snr_sim_MRC);

figure
plot(nRx,snr_sim_log_MRC,'LineWidth',1.5)
xlabel('Number of receive antennas')
ylabel('SNR (dB)')
title('Maximal Ratio Combining (MRC)')
grid on


%% =========================
%  EGC CODE
% =========================

snr_sim_EGC=zeros([length(snr_dB) nRx_max]);

for j=1:nRx_max
    for k=1:length(snr_dB)

        h=randn(j,N)+(randn(j,N)*1i);
        x_kron=kron(ones(nRx(j),1),x);
        c_in=h.*x_kron;

        y=awgn(c_in,snr_dB(k),'measured');

        y_rec = y.*exp(-1i*angle(h));
        y_rec = sum(y_rec,1);

        snr_sim_EGC(k,j) = mean(y_rec.*conj(y_rec))/nRx(j);

    end
end

snr_sim_log_EGC=10*log10(snr_sim_EGC);

figure
plot(nRx,snr_sim_log_EGC,'LineWidth',1.5)
xlabel('Number of receive antennas')
ylabel('SNR (dB)')
title('Equal Gain Combining (EGC)')
grid on


%% =========================
%  SC CODE
% =========================

snr_sim_SC=zeros([length(snr_dB) nRx_max]);

for j=1:nRx_max
    for k=1:length(snr_dB)

        h=randn(j,N)+(randn(j,N)*1i);
        x_kron=kron(ones(nRx(j),1),x);
        c_in=h.*x_kron;

        y=awgn(c_in,snr_dB(k),'measured');

        hPower = h.*conj(h);

        [hMaxVal ,ind] = max(hPower,[],1);

        hMaxValMat = kron(ones(nRx(j),1),hMaxVal);

        y_rec = y(hPower==hMaxValMat);
        hSel = h(hPower==hMaxValMat);

        snr_sim_SC(k,j) = mean(hSel.*conj(hSel));

    end
end

snr_sim_log_SC=10*log10(snr_sim_SC);

figure
plot(nRx,snr_sim_log_SC,'LineWidth',1.5)
xlabel('Number of receive antennas')
ylabel('SNR (dB)')
title('Selection Combining (SC)')
grid on


%% =========================
%  COMPARISON GRAPH
% =========================

MRC_avg = mean(snr_sim_log_MRC,1);
EGC_avg = mean(snr_sim_log_EGC,1);
SC_avg  = mean(snr_sim_log_SC,1);

figure

plot(nRx,MRC_avg,'m','LineWidth',3)
hold on
plot(nRx,EGC_avg,'Color',[1 0.5 0],'LineWidth',3)
plot(nRx,SC_avg,'k','LineWidth',3)

legend('MRC','EGC','SC','Location','southeast')

xlabel('Number of Receiver Antennas')
ylabel('SNR (dB)')
title('MRC - EGC - SC: SNR Versus No. of Receiver Antennas')

grid on


%% =========================
%  TABLE VALUES
% =========================

max_MRC = max(MRC_avg);
mean_MRC = mean(MRC_avg);

max_EGC = max(EGC_avg);
mean_EGC = mean(EGC_avg);

max_SC = max(SC_avg);
mean_SC = mean(SC_avg);

table_data = {
'MRC', sprintf('%.4f',max_MRC), sprintf('%.4f',mean_MRC);
'EGC', sprintf('%.4f',max_EGC), sprintf('%.4f',mean_EGC);
'SC',  sprintf('%.4f',max_SC),  sprintf('%.4f',mean_SC);
};

column_name = {'Technique','Max SNR (dB)','Mean SNR (dB)'};


%% =========================
%  SMALL COMPARISON TABLE
% =========================

uitable('Data',table_data,...
        'ColumnName',column_name,...
        'Position',[45 25 340 90],...
        'FontSize',9,...
        'ColumnWidth',{90 120 120});