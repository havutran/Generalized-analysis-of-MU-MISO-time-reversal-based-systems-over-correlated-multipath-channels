function [R] = CorrMatrix(L,rho)
R=ones(1,L);
for i=1:L
   R(1,i) = R(1,i)*rho^(i-1); 
end

for i=2:L
R(i,:) = circshift (R(i-1,:).',1).';
R(i,1) = R(1,i);
end
end