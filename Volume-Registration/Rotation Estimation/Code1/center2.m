 
 function [P] = center2(P, muP)

n = size(P,1);
%muP = mean(P);
P = P + muP(ones(1,n),:);

return;