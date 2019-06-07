clear;
Mt = 4;
N = 5;
L = 7;
rho1 = 0.5;
rho2 = 0.2;
S=0; 
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

%TR form
for t=0:Mt-1
h_tr(t+1,:) = conj(fliplr(h(2,t*L+1:(t+1)*L)))/sqrt(L);
end

%Equivalent channel
%Psig, Pisi

for t=0:Mt-1
heq (t+1,:) = conv(h_tr(t+1,:), h(1,t*L+1:(t+1)*L) );
end
ht_eq=0;
for t=1:Mt
    ht_eq = ht_eq + heq(t,:);
end

P_heq = abs(ht_eq).^2;
S=S + P_heq;
end
S=S/k;
t=linspace(0,1,2*L-1);
figure(1)
stem (t,S)
F = (sum(S))



%Central tap
a = Mt*(Mt-1)*L*rho2^2 + Mt*(Mt-1)*rho1^2;
%b = 2*(0 + 7*(1*0.6)^2 + 42*(0.6*1)^2)/7;
b = Mt*(L + L*(rho2)^2 + L*(L-1)*(rho2)^2)/L;
c=a+b;

%Total other taps
d=0;
for i = 1: L-1
d = d + 2*i*Mt*1/L + 2*Mt*(Mt-1)*i*rho1^2/L;
end

Ff = c+d






