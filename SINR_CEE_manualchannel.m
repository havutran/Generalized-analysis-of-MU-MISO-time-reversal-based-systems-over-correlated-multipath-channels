clear;
Mt = [2 4 6];
rho1=[0.4];
rho2=[0.5];
N = [2];
L = 6;
PGdb = [0 -3 -10 -18 -26 -32];
psi = linspace(0,2,6);

for k=1:length(Mt)
    k
for i=1:length(psi)
    i
[SINRdb(k,i), SINRdba(k,i),SINR(k,i), SINRa(k,i)] = CEE_manualchannel (Mt(k), N, L, rho1,rho2, psi(1,i),PGdb );
end
end

figure(2); clf;
plot (psi,SINRdb(1,:),'*k',psi,SINRdba(1,:),'-k',psi,SINRdb(2,:),'*k',psi,SINRdba(2,:),'-k',psi,SINRdb(3,:),'*k',psi,SINRdba(3,:),'-k','linewidth',1,'MarkerSize',8)
%grid on;
title('ITU-R Channel Standard - Indoor')
xlabel('\psi')
ylabel('SINR (dBm)')
legend('Simulation','Analytical',1 )

for k=1:length(rho1)
for i=1:length(psi)
reduction(k,i) = 10*log10 (SINR(k,i)/SINR(k,1));
reductiona(k,i) = 10*log10 (SINRa(k,i)/SINRa(k,1));
end
end

figure(3); clf;
plot (psi,reduction(1,:),'*k',psi,reduction(1,:),'-k',psi,reduction(2,:),'*b',psi,reduction(2,:),'-b',psi,reduction(3,:),'*r',psi,reduction(3,:),'-r','linewidth',1,'MarkerSize',8)
%grid on;
title('ITU-R Channel Standard - Indoor, N = 2')
xlabel('\psi')
ylabel('Average Reduction (dB)')
legend('Simulation','Analytical',1 )