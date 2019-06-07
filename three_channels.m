clear;
Mt = 3;
N = 3;
L = 6;
%SNR
SNRdb = linspace(-5,25,31);
SNR = 10.^(SNRdb/10);

%Vehicular
PGdb = [0 -1 -9 -10 -15 -20];
[C(1,:),Ca(1,:),SINRdb(1,:),SINRdba(1,:)] = SINR_manual_channel (Mt,N,L,SNR,PGdb);

%Indoor
PGdb1 = [0 -3 -10 -18 -26 -32];
[Cx1(1,:),Ca1(1,:),SINRdb1(1,:),SINRdba1(1,:)] = SINR_manual_channel (Mt,N,L,SNR,PGdb1);

%FengHan
Ts = 1/(20);
D=1;

%rmdls
PG = 10.^(PGdb/10);
delay = [0 310 710 1090 1730 2510]/1000;
PG1 = 10.^(PGdb1/10);
delay1 = [0 50 110 170 290 310]/1000;

med = delay*PG'/sum(PG);
secondmoment = delay.^2*PG'/sum(PG);
rmsdls = sqrt(secondmoment - med^2);

med1 = delay1*PG1'/sum(PG1);
secondmoment1 = delay1.^2*PG1'/sum(PG1);
rmsdls1 = sqrt(secondmoment1 - med1^2);

%E[Psig]
A = ((Mt+1)*(1 - exp (-(L+1)*Ts/rmsdls)) + (Mt-1)*(exp(-Ts/rmsdls)-exp(-L*Ts/rmsdls)))/(1-exp(-2*Ts/rmsdls));

%E[Pisi]
B = 2*(exp(-Ts/rmsdls)*(1-exp(-(L-2+D)*Ts/rmsdls)))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls)));

%E[Piui]
C = (N-1)*(1+exp(-D*Ts/rmsdls)+exp(-2*L*Ts/rmsdls)-2*exp(-(L+1)*Ts/rmsdls)-2*exp(-(D+L-1)*Ts/rmsdls)+exp(-(D+2*L)*Ts/rmsdls))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls))*(1-exp(-L*Ts/rmsdls)));

%E[Psig]
A1 = ((Mt+1)*(1 - exp (-(L+1)*Ts/rmsdls1)) + (Mt-1)*(exp(-Ts/rmsdls1)-exp(-L*Ts/rmsdls1)))/(1-exp(-2*Ts/rmsdls1));

%E[Pisi]
B1 = 2*(exp(-Ts/rmsdls1)*(1-exp(-(L-2+D)*Ts/rmsdls1)))/((1-exp(-D*Ts/rmsdls1))*(1+exp(-Ts/rmsdls1)));

%E[Piui]
C1 = (N-1)*(1+exp(-D*Ts/rmsdls1)+exp(-2*L*Ts/rmsdls1)-2*exp(-(L+1)*Ts/rmsdls1)-2*exp(-(D+L-1)*Ts/rmsdls1)+exp(-(D+2*L)*Ts/rmsdls1))/((1-exp(-D*Ts/rmsdls1))*(1+exp(-Ts/rmsdls1))*(1-exp(-L*Ts/rmsdls1)));


for i = 1 : length(SNR)
sigma(1,i) = (1/SNR(1,i))*((1-exp(-L*Ts/rmsdls))/(1-exp(-Ts/rmsdls)));
SINR2(1,i) = A/(B+C+sigma(1,i));
SINRdB2(1,i) = 10 *log10(SINR2(1,i));

sigma1(1,i) = (1/SNR(1,i))*((1-exp(-L*Ts/rmsdls1))/(1-exp(-Ts/rmsdls1)));
SINR3(1,i) = A/(B+C+sigma1(1,i));
SINRdB3(1,i) = 10 *log10(SINR3(1,i));
end

figure(3); clf;
plot (SNRdb, SINRdb1,'*k',SNRdb, SINRdba1,'-k',SNRdb, SINRdb,'ok',SNRdb, SINRdba,'-k',SNRdb,SINRdB2, ':k',SNRdb,SINRdB3, '-.k' ,'linewidth',1,'MarkerSize',8)
title('ITU-R Channel Standard')
xlabel('SNR (dBm), M = 3, N = 3')
ylabel('SINR (dBm)')
legend('Simulation - Indoor','Analytical - Indoor','Simulation - Vehicular','Analytical - Vehicular','Han^,s formula [3] - Vehicular','Han^,s formula []3 - Indoor',4 )



%Vehicular
PGdb = [0 -1 -9 -10 -15 -20];
[C5(1,:),Ca5(1,:),SINRdb5(1,:),SINRdba5(1,:)] = SINR_manual_channel (6,N,L,SNR,PGdb);

Mt = 6;
%E[Psig]
A3 = ((Mt+1)*(1 - exp (-(L+1)*Ts/rmsdls)) + (Mt-1)*(exp(-Ts/rmsdls)-exp(-L*Ts/rmsdls)))/(1-exp(-2*Ts/rmsdls));

%E[Pisi]
B3 = 2*(exp(-Ts/rmsdls)*(1-exp(-(L-2+D)*Ts/rmsdls)))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls)));

%E[Piui]
C3 = (N-1)*(1+exp(-D*Ts/rmsdls)+exp(-2*L*Ts/rmsdls)-2*exp(-(L+1)*Ts/rmsdls)-2*exp(-(D+L-1)*Ts/rmsdls)+exp(-(D+2*L)*Ts/rmsdls))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls))*(1-exp(-L*Ts/rmsdls)));


for i = 1 : length(SNR)
sigma2(1,i) = (1/SNR(1,i))*((1-exp(-L*Ts/rmsdls))/(1-exp(-Ts/rmsdls)));
SINR4(1,i) = A3/(B3+C3+sigma2(1,i));
SINRdB4(1,i) = 10 *log10(SINR4(1,i));
end

figure(4); clf;
plot (SNRdb, SINRdb,'*k',SNRdb, SINRdba,'-k',SNRdb,SINRdB2, ':k',SNRdb, SINRdb5,'*k',SNRdb, SINRdba5,'-k',SNRdb,SINRdB4, ':k' ,'linewidth',1,'MarkerSize',8)
title('ITU-R Channel Standard')
xlabel('SNR (dBm), N = 3')
ylabel('SINR (dBm)')
legend('Simulation - Vehicular','Analytical - Vehicular','Han^,s work [5] - Vehicular',4 )



























