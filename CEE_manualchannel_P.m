function [SINRdb, SINRdba,SINR, SINRa] = CEE_manualchannel (Mt, N, L, rho1, rho2, psi, PGdb, SNRdb,K)

%Mt=4;
%rho1=0.8;
%rho2=0.5;
%N = 4;
%L = 6;
%PGdb = [0 -1 -9 -10 -15 -20];
%psi = 1.5;

PG = 10.^(PGdb/10);
sum_p = sum(PG);
x = (1/(1+psi));

S=0;
S1=0;
Sx=0;

%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);

%Simlation
for loop=1:K
%i.i.d channel
hw = sqrt(1/(1+psi))*manual_channel(L,N,Mt,PGdb);
e = sqrt(psi/(1+psi))*manual_channel(L,N,Mt,PGdb);
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
h_tr(t+1,:) = conj(fliplr(h_est(1,t*L+1:(t+1)*L)))/sqrt(Mt*sum_p*x);
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
h_tr1(t+1,:) = conj(fliplr(h_est(2,t*L+1:(t+1)*L)))/sqrt(Mt*sum_p*x);
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

Psig = P_heq(1,L);
Pisi = sum(P_heq) - Psig;
Piui = sum(P_heq1);
SNR = 10.^(SNRdb/10);
%noise = Psig./SNR;
for i=1:length(SNR)
temp(1,i) = SNR(1,i)*Psig/(SNR(1,i)*Pisi + SNR(1,i)*(N-1)*Piui + 1);
end

Sx = Sx + temp;

end
%S=S/k;
%S1=S1/k;
%Psig = S(1,L);
%Pisi = (sum(S) - Psig);
%Piui = sum(S1);

%SNR

%SNR = 10.^(SNRdb/10);
%noise = Psig./SNR;
%for i=1:length(SNR)
%SINR(1,i) = Psig/(Pisi + (N-1)*Piui + noise(1,i));
%end
SINR = Sx/loop;
SINRdb = 10.*log10 (Sx/loop);

%Analytical


%Pisi
S1a=0;
w=0;
xi_tong1=0;
for k = 1: L-1
    for l=1:k
       w = w + PG(k+1-l)*PG(L+1-l); 
    end
    for l=1:k
       xi_tong1 = xi_tong1 + PG(k+1-l)*PG(L+1-l)*x*rho1^2; 
    end
end

S1a = S1a + 2*w/sum_p + 2*(Mt-1)*xi_tong1*(1+psi)/sum_p;
Pisia = S1a;

%Psig
xi_tong2 = sum(PG.^2)*psi*x*rho1^2;
%Psiga = ((L + L^2)/L + (Mt-1)*L+ (Mt-1)*rho1^2*sqrt(1));
Psiga = (sum(PG.^2) + (sum_p)^2*x)/(sum_p) + (Mt-1)*x*((sum_p)^2 + sum(PG.^2)*rho1^2 + (1+psi)*xi_tong2)/sum_p;

%Piui
%Central tap
w1=0;
%k = L;
xi_tong3 = sum(PG.^2)*(rho2^2*x^2 + rho1^2*x);
xi_tong4 = 0;

for l1=1:L
    for l2=1:L
        if l2 ~= l1
            xi_tong4 = xi_tong4 + PG(l1)*PG(l2)*x^2*rho2^2;
        end
    end
end

P_central = ((1+rho2^2*x)*sum(PG.^2) + xi_tong4*(1+psi))/sum_p + (1+psi)*(Mt-1)*(xi_tong4 + xi_tong3)/sum_p;

%for l=1:k
%   w1 = w1 + PG(k+1-l)*PG(L+1-l); 
%end


%Total other taps
%d=0;
%for i = 1: L-1
%d = d + 2*i*1/L + 2*(Mt-1)*i*rho1^2/L;
%end

Piuia = Pisia + P_central;

for i=1:length(SNR)
SINRa(1,i) = SNR(1,i)*Psiga/(SNR(1,i)*Pisia + SNR(1,i)*(N-1)*Piuia + 1);
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




























