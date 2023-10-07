 
 function [P,muP] = center(P)

n = size(P,1);
 muP = mean(P);
%muP=[P(7,1),P(7,2),P(7,3)]
P = P - muP(ones(1,n),:);

return;