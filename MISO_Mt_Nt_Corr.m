clear;
for k=1:length(Mt)
N = 4;
L = 7;
rho1 = linspace(0,0.9,10);
rho2 = linspace(0,0.9,10);
for i=1:length(rho1)
[SINRdb(k,i), SINRdba(k,i),SINR(k,i), SINRa(k,i)] = SINR_TR_MU_MISO (Mt(1,k), N, L, rho1(1,i), rho2(1,1));
end

end
























