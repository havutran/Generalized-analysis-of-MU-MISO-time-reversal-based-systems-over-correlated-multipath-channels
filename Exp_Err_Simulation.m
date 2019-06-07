clear;
%Simulation

Mt=4;
rho1=0;
rho2=0;
N = 2;
L = 6;
PGdb = [0 -1 -9 -10 -15 -20];
psi = 1.5;

x = (1/(1+psi));

S=0;
S1=0;

%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);

hw = manual_channel(L,N,Mt,PGdb);
sum_p1 = norm(hw(1,:))^2;
sum_p2 = norm(hw(2,:))^2;

hwdb = 10.*log10 (abs(hw).^2);
for k=1:10000
%i.i.d channel

%e = sqrt(psi/(1+psi))*manual_channel(L,N,Mt,PGdb);

for n=1:N
    for t=0:Mt-1
     e(n,t*L+1 : (t+1)*L) =  sqrt(1/(1+psi))*manual_channel(L,1,1,hwdb(n,t*L+1 : (t+1)*L));
    end
end

h_p = e;
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
h_tr(t+1,:) = conj(fliplr(h_est(1,t*L+1:(t+1)*L)))/sqrt(sum_p1);
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
h_tr1(t+1,:) = conj(fliplr(h_est(2,t*L+1:(t+1)*L)))/sqrt(sum_p2);
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
SNRdb=5;
SNR = 10.^(SNRdb/10);
noise = Psig./SNR;
for i=1:length(SNR)
SINR(1,i) = Psig/(Pisi + (N-1)*Piui + noise(1,i));
end
SINRdb = 10.*log10 (SINR);