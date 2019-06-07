clear;

Mt = [2];
rho1=[0];
rho2=[0];
N = [2];
L = 110;
%B = 500MHz
Ts = 10^-9;
rmsdls = 100*Ts;
psi = 0;
for s = 1 : L
    PGdb (1,s) = 1 * exp (-Ts * (s-1)/rmsdls);      
end

%SNR
SNRdb = linspace(-5,25,16);
SNR = 10.^(SNRdb/10);

%Our analysis
[C1(1,:),Ca1(1,:),SINRdb1(1,:),SINRdba1(1,:)] = SINR_manual_channel (Mt,N,L,SNR,PGdb);

%Feng Han
[SINR3, SINRdB3] = Han_analysis (L,N,Mt,SNR,Ts, rmsdls, PGdb);

Mt = 4;
[SINR4, SINRdB4] = Han_analysis (L,N,Mt,SNR,Ts, rmsdls, PGdb);

%Our analysis
[C2(1,:),Ca2(1,:),SINRdb2(1,:),SINRdba2(1,:)] = SINR_manual_channel (Mt,N,L,SNR,PGdb);

Mt=6;
[C5(1,:),Ca5(1,:),SINRdb5(1,:),SINRdba5(1,:)] = SINR_manual_channel (Mt,N,L,SNR,PGdb);
[SINR6, SINRdB6] = Han_analysis (L,N,Mt,SNR,Ts, rmsdls, PGdb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Corr + CEE
Mt = 2;
[SINRdb7, SINRdba7,SINR7, SINRa7] = CEE_manualchannel (Mt, N, L, 0.4, 0.5, 0.5, PGdb, SNRdb,1000);
Mt = 4;
[SINRdb8, SINRdba8,SINR8, SINRa8] = CEE_manualchannel (Mt, N, L, 0.4, 0.5, 0.5, PGdb, SNRdb,1000);
Mt = 6;
[SINRdb9, SINRdba9,SINR9, SINRa9] = CEE_manualchannel (Mt, N, L, 0.4, 0.5, 0.5, PGdb, SNRdb,1000);


figure(4); clf;
plot (SNRdb, SINRdb1,'*r',SNRdb, SINRdba1,'-b',SNRdb,SINRdB3, 'ok',SNRdb, SINRdb2,'*r',SNRdb, SINRdba2,'-b',SNRdb,SINRdB4, 'ok',SNRdb, SINRdb5,'*r',SNRdb, SINRdba5,'-b',SNRdb,SINRdB6, 'ok' ,SNRdb, SINRdb7,'*r',SNRdb, SINRdba7,'-b',SNRdb, SINRdb8,'*r',SNRdb, SINRdba8,'-b',SNRdb, SINRdb9,'*r',SNRdb, SINRdba9,'-b','linewidth',1,'MarkerSize',8)
title('Wideband Channel')
xlabel('SNR (dBm), N = 3')
ylabel('SINR (dBm)')
legend('Simulation','Analytical - Our fomulation','Analytical - Han^,s work [5]',4 )





