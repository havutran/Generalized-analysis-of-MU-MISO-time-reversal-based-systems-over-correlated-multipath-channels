clear;
Mt = [6];
rho1=[0 0.5 0.9];
rho2=[0 0.5 0.9];
N = [2];
L = 110;
%B = 500MHz
Ts = 10^-9;
rmsdls = 100*Ts;

for s = 1 : L
    PGdb (1,s) = 1 * exp (-Ts * (s-1)/rmsdls);      
end
%PGdb = [0 -3 -10 -18 -26 -32];
psi = linspace(0,2,6);

for k=1:length(rho1)
    k
for i=1:length(psi)
    i
[SINRdb(k,i), SINRdba(k,i),SINR(k,i), SINRa(k,i)] = CEE_manualchannel (Mt, N, L, rho1(k),rho2(k), psi(1,i),PGdb );
end
end

figure(1); clf;
plot (psi(1,1), SINRdba (1,1),'or' ,psi,SINRdb(1,:),'*k',psi,SINRdba(1,:),'-k',psi,SINRdb(2,:),'*k',psi,SINRdba(2,:),'-k',psi,SINRdb(3,:),'*k',psi,SINRdba(3,:),'-k','linewidth',1,'MarkerSize',8)
%grid on;
title('Wideband Channel, M = 6, N = 2')
xlabel('\psi (error)')
ylabel('SINR (dB)')
legend('Analytical - Formulas of Han et al. []','Simulation','Analytical',1 )

