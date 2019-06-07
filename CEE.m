function [SINRdb, SINRdba,SINR, SINRa] = CEE (Mt, N, L, rho1,rho2, psi )
S=0;
S1=0;

%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);

%Simlation
for k=1:10000
%i.i.d channel
hw = sqrt(1/(1+psi))*(randn(N,Mt*L) + 1j*randn(N,Mt*L))/sqrt(2);
e = sqrt(psi/(1+psi))*(randn(N,Mt*L) + 1j*randn(N,Mt*L))/sqrt(2);
h_p = hw + e;
for t=0:Mt-1
hw1 ( t*N+1 : (t+1)*N  , 1:L  ) = hw ( 1:N    ,  t*L+1 : (t+1)*L  );
h_p1 ( t*N+1 : (t+1)*N  , 1:L  ) = h_p ( 1:N    ,  t*L+1 : (t+1)*L  );
end

%Correlated channel
corr = kron((R_t.')^(1/2),R_u^(1/2));
h_est1 = corr*hw1;
h_pcorr1 = corr*h_p1;
%h = reshape (h.', Mt*L, N).';
for t=0:Mt-1
h_est ( 1:N    ,  t*L+1 : (t+1)*L  )= h_est1 ( t*N+1 : (t+1)*N  , 1:L  );
h_pcorr ( 1:N    ,  t*L+1 : (t+1)*L  )= h_pcorr1 ( t*N+1 : (t+1)*N  , 1:L  );
end

%Perfect channel
%h_p = h_est + e;

%Equivalent channel
%Psig, Pisi
%TR form
for t=0:Mt-1
h_tr(t+1,:) = conj(fliplr(h_est(1,t*L+1:(t+1)*L)))/sqrt(Mt*L/(1+psi));
end

for t=0:Mt-1
heq (t+1,:) = conv(h_tr(t+1,:), h_pcorr(1,t*L+1:(t+1)*L) );
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
h_tr1(t+1,:) = conj(fliplr(h_est(2,t*L+1:(t+1)*L)))/sqrt(Mt*L/(1+psi));
end
for t=0:Mt-1
heq1 (t+1,:) = conv(h_tr1(t+1,:), h_pcorr(1,t*L+1:(t+1)*L) );
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
x = (1/(1+psi));
%Psig Pisi
S1a=0;
for i = 1: L-1
    S1a = S1a + 2*i*( x^2 + psi*x^2 )/(L*x)+ 2*(Mt-1)*rho1^2*i*x/(L*x)  ;
end
Pisia = S1a;

Psiga = ((L*x^2 + (L*x)^2 + L*x^2*psi)/(L*x) + (Mt-1)*(L*x)+ (Mt-1)*L*x^2*(1+(psi))*rho1^2/(L*x));

% 21.4578
%Psiga = ((L + L^2)/L + (Mt-1)*L+ (Mt-1)*rho1^2*sqrt(1));

%Piui
%Central tap

a =  (Mt-1)*(rho1^2*(x)^2*psi*L + rho2^2*x^2*L^2 + L*rho1^2*x^2  )/(L*x);
%b = 2*(0 + 7*(1*0.6)^2 + 42*(0.6*1)^2)/7;
b = (L*x + L^2*x^2*rho2^2 )/(L*x);
c=a+b;


%Total other taps
d=0;
for i = 1: L-1
d = d + 2*i*(1+psi)*x^2/(L*x) + 2*(Mt-1)*i*(rho1^2*x^2+ (psi)*(x^2)*rho1^2 )/(L*x);
end

Piuia = c+d;
for i=1:length(SNR)
SINRa(1,i) = Psiga/(Pisia + (N-1)*Piuia + noise(1,i));
end
SINRdba = 10.*log10 (SINRa);

end

%P=0;P1=0;P2=0;
%for i = 1:10000
%h_est = sqrt(1/1.2)*(randn(1,6) + 1j*randn(1,6))/sqrt(2);
%e = sqrt(0.2/1.2)*(randn(1,6) + 1j*randn(1,6))/sqrt(2);
%h = h_est+e;
%P = P + norm(h)^2;
%P1 = P1 + norm(h_est)^2;
%P2 = P2 + norm(e)^2;
%end
%P=P/i
%P1=P1/i
%P2=P2/i
%P1+P2




























