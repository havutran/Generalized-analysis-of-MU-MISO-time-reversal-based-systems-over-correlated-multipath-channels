function [SINRdb, SINRdba,SINR, SINRa] = SINR_TR_MU_MISO (Mt, N, L, rho1, rho2)
%Simulation

S=0;S1=0; 
%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);
%Simlation
for k=1:10000
%i.i.d channel
hw = (randn(N,Mt*L) + 1j*randn(N,Mt*L))/sqrt(2);

for t=0:Mt-1
hw1 ( t*N+1 : (t+1)*N  , 1:L  ) = hw ( 1:N    ,  t*L+1 : (t+1)*L  );
end

%Correlated channel
corr = kron((R_t.')^(1/2),R_u^(1/2));
h = corr*hw1;
%h = reshape (h.', Mt*L, N).';
for t=0:Mt-1
h ( 1:N    ,  t*L+1 : (t+1)*L  )= h ( t*N+1 : (t+1)*N  , 1:L  );
end

%Equivalent channel
%Psig, Pisi
%TR form
for t=0:Mt-1
h_tr(t+1,:) = conj(fliplr(h(1,t*L+1:(t+1)*L)))/sqrt(Mt*L);
end
for t=0:Mt-1
heq (t+1,:) = conv(h_tr(t+1,:), h(1,t*L+1:(t+1)*L) );
end
ht_eq=0;
for t=1:Mt
    ht_eq = ht_eq + heq(t,:);
end
P_heq = abs(ht_eq).^2;
S=S + P_heq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Piui
%TR form
for t=0:Mt-1
h_tr1(t+1,:) = conj(fliplr(h(2,t*L+1:(t+1)*L)))/sqrt(Mt*L);
end
for t=0:Mt-1
heq1 (t+1,:) = conv(h_tr1(t+1,:), h(1,t*L+1:(t+1)*L) );
end
ht_eq1=0;
for t=1:Mt
    ht_eq1 = ht_eq1 + heq1(t,:);
end

P_heq1 = abs(ht_eq1).^2;
S1=S1 + P_heq1;

end
S=S/k;
S1=S1/k;
Psig = S(1,L);
Pisi = (sum(S) - Psig);
Piui = sum(S1);

%SNR
SNRdb = 5;
SNR = 10.^(SNRdb/10);
noise = Psig./SNR;
for i=1:length(SNR)
SINR(1,i) = Psig/(Pisi + (N-1)*Piui + noise(1,i));
end
SINRdb = 10.*log10 (SINR);

%Analytical
%Psig Pisi
S1a=0;
for i = 1: L-1
    S1a = S1a + 2*i/L  + 2*(Mt-1)*(rho1^2*i)/(L)  ;
end
Pisia = S1a;

Psiga = ((L + L^2)/L + (Mt-1)*L+ (Mt-1)*rho1^2*sqrt(1));

%Piui
%Central tap
a = (Mt-1)*L*rho2^2 + (Mt-1)*rho1^2;
%b = 2*(0 + 7*(1*0.6)^2 + 42*(0.6*1)^2)/7;
b = (L + L*(rho2)^2 + L*(L-1)*(rho2)^2)/L;
c=a+b;

%Total other taps
d=0;
for i = 1: L-1
d = d + 2*i*1/L + 2*(Mt-1)*i*rho1^2/L;
end

Piuia = c+d;

for i=1:length(SNR)
SINRa(1,i) = Psiga/(Pisia + (N-1)*Piuia + noise(1,i));
end
SINRdba = 10.*log10 (SINRa);

end











