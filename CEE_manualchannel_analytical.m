function [SINRax,SINRdba,SINRa, Psiga, Pisia, P_central, noise] = CEE_manualchannel_analytical (Mt, N, L, rho1, rho2, psi, PGdb, SNRdb)

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

%Correlation matrix
R_t = CorrMatrix_interclass(Mt,rho1);
R_u = CorrMatrix_interclass(N,rho2);

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


Piuia = Pisia + P_central;


SNR = 10.^(SNRdb/10);

for i=1:length(SNR)
SINRa(1,i) = SNR(i)*Psiga/(SNR(i)*Pisia + SNR(i)*(N-1)*Piuia + 1);
end
noise=1;
SINRax = Psiga/(Pisia + (N-1)*Piuia + 1);
SINRdba = 10.*log10 (SINRa);


end






























