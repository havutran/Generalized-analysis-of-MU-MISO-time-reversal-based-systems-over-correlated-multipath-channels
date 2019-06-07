clear;
%For transmitter side
Mt = [2 4 6];
for k=1:length(Mt)
N = 4;
L = 7;
rho1 = linspace(0,1,11);
rho2 = linspace(0,1,11);
for i=1:length(rho1)
[SINRdb(k,i), SINRdba(k,i),SINR(k,i), SINRa(k,i)] = SINR_TR_MU_MISO (Mt(1,k), N, L, rho1(1,i), rho2(1,3));
end

end

for k = 1:length(Mt)
for i = 1: length(rho1)
Gain(k,i) = 10*log10 (SINR(k,i)/SINR(k,1));
end
end

for k = 1:length(Mt)
for i = 1: length(rho1)
capacity(k,i) = log2 (1 + SINR(k,i));
capacitya(k,i) = log2 (1 + SINRa(k,i));
end
end


%For receiver side
N = [2 4 6];
for k=1:length(N)
Mt = 4;
L = 117;
rho1 = linspace(0,1,11);
rho2 = linspace(0,1,11);
for i=1:length(rho2)
[SINRdb_r(k,i), SINRdb_ra(k,i),SINR_r(k,i), SINR_ra(k,i)] = SINR_TR_MU_MISO (Mt, N(1,k), L, rho1(1,3), rho2(1,i));
end

end

for k = 1:length(N)
for i = 1: length(rho2)
Gain_r(k,i) = 10*log10 (SINR_ra(k,i)/SINR_ra(k,1));
end
end

for k = 1:length(N)
for i = 1: length(rho2)
capacity_r(k,i) = log2 (1 + SINR_r(k,i));
capacity_ra(k,i) =log2 (1 + SINR_ra(k,i));
end
end

figure(1); clf;
plot (rho1, SINRdb(1,:),'*k',rho1, SINRdba(1,:),'-k',rho1, SINRdb(2,:),'ok',rho1, SINRdba(2,:),'--k',rho1, SINRdb(3,:),'^k',rho1, SINRdba(3,:),':k','linewidth',1,'MarkerSize',8)
%grid on;
title('')
xlabel('\rho_T')
ylabel('SINR (dB)')
legend('Simulation','Analytical',1 )


figure(2); clf;
plot (rho1, Gain(1,:),'-k',rho1, Gain(2,:),'--k',rho1, Gain(3,:),':k',rho1, Gain_r(1,:),'-*k',rho1, Gain_r(2,:),'--*k',rho1, Gain_r(3,:),':*k','linewidth',1,'MarkerSize',8)
%grid on;
title('')
xlabel('\rho_T')
ylabel('Gain (dB)')
legend('M_T = 2','M_T = 4','M_T = 6',1 )


figure(3); clf;
plot (rho1, capacity(1,:),'*k',rho1, capacitya(1,:),'-k',rho1, capacity(2,:),'ok',rho1, capacitya(2,:),'--k',rho1, capacity(3,:),'^k',rho1, capacitya(3,:),':k','linewidth',1,'MarkerSize',8)
%grid on;
title('')
xlabel('\rho_T')
ylabel('Capacity')
legend('Simulation','Analytical',1 )

figure(4); clf;
plot (rho1, capacity_r(1,:),'*k',rho1, capacity_ra(1,:),'-k',rho1, capacity_r(2,:),'ok',rho1, capacity_ra(2,:),'--k',rho1, capacity_r(3,:),'^k',rho1, capacity_ra(3,:),':k','linewidth',1,'MarkerSize',8)
%grid on;
title('')
xlabel('\rho_U')
ylabel('Capacity')
legend('Simulation','Analytical',1 )

