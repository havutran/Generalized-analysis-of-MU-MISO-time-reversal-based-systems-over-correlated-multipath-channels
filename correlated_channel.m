clear;
m=2;
L=6;
rho = 0.3;
S1=0; S2=0;
%Correlation matrix
R_t = CorrMatrix(m,rho);
R_t1 = CorrMatrix(6,rho);
for k=1:100000
hw(1,:) = (randn(1,L) + 1j*randn(1,L))/sqrt(2);
hw(2,:) = (randn(1,L) + 1j*randn(1,L))/sqrt(2);

%Correlated channel
h = R_t^(1/2)*hw*R_t1^(1/2);
S1 = S1 + norm(h(1,:))^2;
S2 = S2 + norm(h(2,:))^2;
end

N1 = S1/k
N2 = S2/k