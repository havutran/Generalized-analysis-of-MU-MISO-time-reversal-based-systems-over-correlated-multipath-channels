

Mt=4;
rho1=0;
rho2=0;
N = 2;
L = 6;
psi = 1.5;


sum_p1 = norm(hw(1,:))^2;
sum_p2 = norm(hw(2,:))^2;

x = (1/(1+psi));
%Analytical


%Pisi
xi_tong1 = 0;
for k =1:L-1
for t=0:Mt-1
for l1 = 1:k
            xi_tong1 = xi_tong1 + 2*abs(hw(1,t*L+1+k-l1))^2*abs(hw(1,t*L+1+L-l1))^2*x/sum_p1;
end
end
end
Pisia = xi_tong1;

%Psig

Psiga = (x*sum(abs(hw(1,:)).^4) )/(sum_p1) ;

%Piui
%Central tap
xi_tong3=0;
for l=1:Mt*L
xi_tong3 = xi_tong3 + x*abs(hw(1,l))^2*abs(hw(2,l))^2/sum_p2;
end
Piui_central = xi_tong3;

%Othertap
xi_tong4 = 0;
for k =1:L-1
for t=0:Mt-1
for l1 = 1:k
            xi_tong4 = xi_tong4 + abs(hw(1,t*L+1+k-l1))^2*abs(hw(2,t*L+1+L-l1))^2*x/sum_p2;
end
end
end

xi_tong6 = 0;
for k =1:L-1
for t=0:Mt-1
for l1 = 1:k
            xi_tong6 = xi_tong6 + abs(hw(1,t*L+ L - (1+k-l1) + 1))^2*abs(hw(2,t*L+ L - (1+L-l1) + 1))^2*x/sum_p2;
end
end
end

Piui_other = xi_tong4+xi_tong6;


Piuia = Piui_central + Piui_other;

for i=1:length(SNR)
SINRa(1,i) = Psiga/(Pisia + (N-1)*Piuia + noise(1,i));
end

SINRdba = 10.*log10 (SINRa);

