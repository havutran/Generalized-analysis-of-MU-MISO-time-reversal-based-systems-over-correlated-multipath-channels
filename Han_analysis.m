function [SINR, SINRdB] = Han_analysis (L,N,Mt,SNR,Ts, rmsdls, PGdb)

D=1;
%E[Psig]
A = ((Mt+1)*(1 - exp (-(L+1)*Ts/rmsdls)) + (Mt-1)*(exp(-Ts/rmsdls)-exp(-L*Ts/rmsdls)))/(1-exp(-2*Ts/rmsdls));

%E[Pisi]
B = 2*(exp(-Ts/rmsdls)*(1-exp(-(L-2+D)*Ts/rmsdls)))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls)));

%E[Piui]
C = (N-1)*(1+exp(-D*Ts/rmsdls)+exp(-2*L*Ts/rmsdls)-2*exp(-(L+1)*Ts/rmsdls)-2*exp(-(D+L-1)*Ts/rmsdls)+exp(-(D+2*L)*Ts/rmsdls))/((1-exp(-D*Ts/rmsdls))*(1+exp(-Ts/rmsdls))*(1-exp(-L*Ts/rmsdls)));


for i = 1 : length(SNR)
sigma(1,i) = (1/SNR(1,i))*Mt*((1-exp(-L*Ts/rmsdls))/(1-exp(-Ts/rmsdls)));
SINR(1,i) = A/(B+C+sigma(1,i));
SINRdB(1,i) = 10 *log10(SINR(1,i));
end

end