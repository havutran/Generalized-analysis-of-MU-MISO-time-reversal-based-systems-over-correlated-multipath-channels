function [R] = CorrMatrix_interclass (L,rho)

R = rho*ones(L,L);

for i=1:L
    R(i,i) = 1;
end
end