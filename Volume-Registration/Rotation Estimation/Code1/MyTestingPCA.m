clear 
clc
close all

a = 10*randn(1000,1);
b = 2*randn(1000,1);
figure, plot (a,b, '.')

R = [ cosd(75) -sind(75);sind(75) cosd(75)];


G = R*[a,b]';

size(G)

x = G(1,:);
y = G(2,:);

hold on
plot(x,y, '.g')
hold off

R

CovMat=cov([x',y']);
[V,D]=eigs(CovMat,2);


EstimateR = fliplr(V)