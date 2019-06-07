clear;
Mt = [6];
rho1=[0 0.4 0.6];
rho2=[0 0.4 0.6];
N = [3];
L = 18;
PGdb = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18 -22.4 -26.7];
%psi = [0 1 2];
%rho = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
psi = linspace(0,3,6);

Ts = 10^-9;
rmsdls = 100*Ts;
for s = 1 : 110
    PGdb1 (1,s) = 1 * exp (-Ts * (s-1)/rmsdls);      
end
%PGdb1 = [0 -0.9 -1.7 -2.6 -3.5 -4.3 -5.2 -6.1 -6.9 -7.8 -4.7 -7.3 -9.9 -12.5 -13.7 -18 -22.4 -26.7];

for k=1:length(rho1)
    k
for i=1:length(psi)
    i
[SINRdb(k,i), SINRdba(k,i),SINR(k,i), SINRa(k,i)] = CEE_manualchannel (Mt, N, 18, rho1(k),rho2(k), psi(i),PGdb, 1, 10000 );
[SINRdb1(k,i), SINRdba1(k,i),SINR1(k,i), SINRa1(k,i)] = CEE_manualchannel (Mt, N,110, rho1(k),rho2(k), psi(i),PGdb1,1,1000);
end
end

figure(2); clf;
plot (psi,SINRdb(1,:),'*r',psi,SINRdba(1,:),'-k',psi,SINRdb1(1,:),'*b',psi,SINRdba1(1,:),'-r',psi,SINRdb(2,:),'*r',psi,SINRdba(2,:),'-k',psi,SINRdb(3,:),'*r',psi,SINRdba(3,:),'-k',psi,SINRdb1(2,:),'*b',psi,SINRdba1(2,:),'-r',psi,SINRdb1(3,:),'*b',psi,SINRdba1(3,:),'-r','linewidth',1,'MarkerSize',8)
%plot (rho,SINRdb(1,:),'*r',rho,SINRdba(1,:),'-k',rho,SINRdb1(1,:),'or',rho,SINRdba1(1,:),'-r',rho,SINRdb(2,:),'*r',rho,SINRdba(2,:),'-k',rho,SINRdb(3,:),'*r',rho,SINRdba(3,:),'-k',rho,SINRdb1(1,:),'*r',rho,SINRdba1(1,:),'-r',rho,SINRdb1(2,:),'*r',rho,SINRdba1(2,:),'-r',rho,SINRdb1(3,:),'*r',rho,SINRdba1(3,:),'-r','linewidth',1,'MarkerSize',8)
%grid on;
%plot (psi,SINRdb1(1,:),'*b',psi,SINRdba1(1,:),'-r',psi,SINRdb1(2,:),'*b',psi,SINRdba1(2,:),'-r','linewidth',1,'MarkerSize',8)
title('UWB and WLAN channels')
xlabel('\psi')
ylabel('Average SINR (dB)')
legend('Simulation - WLAN','Analytical','Simulation - UWB',1 )


%[SINRaxx,SINRdbax,SINRax, Psiga, Pisia, P_central, noise] = CEE_manualchannel_analytical (Mt, N, L, 0, 0, 0, PGdb, 5)
