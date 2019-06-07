clear;
Mt=4;
rho1=0.4;
rho2=0.6;
N = 4;
L = 7;
psi = linspace(0,2,5);

for i=1:length(psi)
[SINRdb(1,i), SINRdba(1,i),SINR(1,i), SINRa(1,i)] = CEE (Mt, N, L, rho1,rho2, psi(1,i) );
end

figure(1); clf;
plot (psi,SINRdb,'*k',psi,SINRdba,'-k','linewidth',1,'MarkerSize',8)
%grid on;
title('\rho_T=0.4,\rho_U=0.6,M_T=4,N=4,L=7')
xlabel('\psi')
ylabel('SINR (dB)')
legend('Simulation','Analytical',1 )

