clear
close all
clc

%% Random outliers on a perfect grid
m = 100; % num points
n = 65;
x=round(rand(1,m)*n);
y=round(rand(1,m)*n);
z=round(rand(1,m)*n);
figure,
scatter3(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis equal