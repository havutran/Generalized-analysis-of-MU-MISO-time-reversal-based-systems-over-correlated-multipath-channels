function [hd] = manual_channel (L,N,Mt,PGdb)

%Number of taps
%L = 6;
%N = 4;
%Mt = 3;

%Average path gain
%PGdb = [0 -3 -10 -18 -26 -32];
PG = 10.^(PGdb/10);
%S=0;
%for k=1:10000
%Random Channel
h = (randn(N,Mt*L) + 1j*randn(N,Mt*L))/sqrt(2);

%Desired Channel
for t=1:Mt
for l = (t-1)*L+1:(t)*L
    for n=1:N
       hd(n,l) = sqrt(PG(l-(t-1)*L))*h(n,l);
    end
end
end
%H = abs(hd).^2;
%S = S + H;
%end



