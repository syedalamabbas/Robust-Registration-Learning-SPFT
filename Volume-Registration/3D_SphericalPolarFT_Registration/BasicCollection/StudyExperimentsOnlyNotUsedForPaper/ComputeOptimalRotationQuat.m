function [ testObjectVertices,M ] = ComputeOptimalRotationQuat( testObjectVertices,referenceObjectVertices )
%COMPUTEOPTIMALROTATIONQUAT This function is a general procedure to compute
%optimal rotation between two object surfaces so long as all the points are
%available -- We only bother about rotation since we are using Fourier
% magnitude data
%================================================================
% This is my attempt to derive the solution proposed in 2013 Transactions
% in Image Processing by Salah and Mahoor

fprintf(' running optimal rotation algorithm \n');

M = eye(4);

% referenceObjectVertices = circshift(referenceObjectVertices,2000);      % Vertices position dependent rotation estimation

% center P and X first
[testObjectVertices,cmP] = center(testObjectVertices);
[referenceObjectVertices,cmX] = center(referenceObjectVertices);

Xc = referenceObjectVertices;

%muP   = mean(testObjectVertices);
%muXc  = mean(Xc);
N     = length(testObjectVertices);
sig   = zeros(3);
for k = 1 : N
    sig = sig + (testObjectVertices(k,:)' * Xc(k,:));
end
sig   = sig/N;
%sig   = sig - (muP' * muXc); % cross-covariance matrix


A     = sig - sig';
delta = [A(2,3) ; A(3,1) ; A(1,2)];
Q     = [trace(sig)  delta' ; delta  sig+sig'-trace(sig)*eye(3)];
[V,D] = eig(Q);
[val,ind] = max(diag(D));
qR    = V(:,ind); % optimal rotation
q0    = qR(1);
q1    = qR(2);
q2    = qR(3);
q3    = qR(4);
R     = [q0^2+q1^2-q2^2-q3^2  2*(q1*q2 - q0*q3)  2*(q1*q3 + q0*q2) ; ...
         2*(q1*q2 + q0*q3)  q0^2+q2^2-q1^2-q3^2  2*(q2*q3 - q0*q1) ; ...
         2*(q1*q3 - q0*q2) 2*(q2*q3 + q0*q1)  q0^2+q3^2-q1^2-q2^2];

Pnew  = (R*testObjectVertices')';
testObjectVertices = Pnew;

M(1:3,1:3) = real(R);

return;

%
% center point set
%
function [P,muP] = center(P)
n = size(P,1);
muP = mean(P);
P = P - muP(ones(1,n),:);
return;



