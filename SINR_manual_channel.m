function [Ct,Ca,SINRdb,SINRdba] = SINR_manual_channel (Mt,N,L,SNR,PGdb)
%Simulation
PG = 10.^(PGdb/10);
sum_p = sum(PG);
rho1 = 0;
rho2 = 0;
S=0;S1=0; Ct=0;
%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);
%Simlation
for k=1:10000
%i.i.d channel
hw = manual_channel(L,N,Mt,PGdb);

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
h_tr(t+1,:) = conj(fliplr(h(1,t*L+1:(t+1)*L)))/sqrt(Mt*sum_p);
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
h_tr1(t+1,:) = conj(fliplr(h(2,t*L+1:(t+1)*L)))/sqrt(Mt*sum_p);
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

noise1 = P_heq(1,L)./SNR;
for i=1:length(SNR)
C(1,i) = P_heq(1,L)/(  (sum(P_heq )- P_heq(1,L)  ) + (N-1)*sum(P_heq1) + noise1(1,i));
end


Ct = Ct + C;

end

Ct = Ct/k;

S=S/k;
S1=S1/k;
Psig = S(1,L);
Pisi = (sum(S) - Psig);
Piui = sum(S1);


noise = Psig./SNR;
for i=1:length(SNR)
SINR(1,i) = Psig/(Pisi + (N-1)*Piui + noise(1,i));
end
SINRdb = 10.*log10 (Ct);

%Analytical
%Psig Pisi
S1a=0;
w=0;
for k = 1: L-1
    for l=1:k
       w = w + PG(k+1-l)*PG(L+1-l); 
    end
end
S1a = S1a + 2*w/sum_p;
Pisia = S1a;

%Psiga = ((L + L^2)/L + (Mt-1)*L+ (Mt-1)*rho1^2*sqrt(1));
Psiga = (sum(PG.^2) + sum_p^2)/sum_p + (Mt-1)*sum_p;

%Piui
%Central tap
w1=0;
k = L;
    for l=1:k
       w1 = w1 + PG(k+1-l)*PG(L+1-l); 
    end


%Total other taps
d=0;
for i = 1: L-1
d = d + 2*i*1/L + 2*(Mt-1)*i*rho1^2/L;
end

Piuia = (2*w+w1)/sum_p;

for i=1:length(SNR)
SINRa(1,i) = Psiga/(Pisia + (N-1)*Piuia + noise(1,i));
end

SINRdba = 10.*log10 (SINRa);


for i=1:length(SNR)
Ca(i) = log2 (1+SINRa(i));
end

end


