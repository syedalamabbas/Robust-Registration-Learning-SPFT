

close all

clear
clc

%% Load and format data
% Let f be the function whose samples are in randomS2sigA bw8.dat, and h be the function whose samples
% are in randomS2sig bw8.dat.We wish to determine how to rotate h so that the correlation is maximized.We
% can think of this graphically: how do we rotate h so that its graph matches that of f’s

Image = load('randomS2sigA_bw8.dat');     % <-------Signal f
Image2 = load('randomS2sig_bw8.dat');     % <-------Pattern h

EditedImage = Image((2*(0:255)+1));
EditedImage2 = Image2((2*(0:255)+1));

N = length(EditedImage);
Final2DImage = reshape(EditedImage, [sqrt(N) sqrt(N)]) ;
Final2DImage2 = reshape(EditedImage2, [sqrt(N) sqrt(N)]);

%% Display

% figure,
% subplot(1,2,1)
% imagesc(Final2DImage)
% xlabel('x')
% ylabel('y')
% subplot(1,2,2)
% imagesc(Final2DImage2)
% xlabel('x')
% ylabel('y')

%% Run the function MEX
addpath(genpath('..\MEXProject'))
fprintf('\n\n\nThis is the old mex project\n\n\n')
[alpha, beta, gamma] = MEXProject(EditedImage, zeros(size(EditedImage)), EditedImage2, zeros(size(EditedImage2)),8,8,7);


%% New SOFT MEX project for arbitrary input sizes
addpath(genpath('..\..\..\MEX_SOFT_Project'))
fprintf('\n\n\nThis is the new mex project\n\n\n')
[alpha, beta, gamma] = MEX_SOFT_Project(0,EditedImage, zeros(size(EditedImage)), EditedImage2, zeros(size(EditedImage2)),8,8,7);

%% Rotating a function and estimating its rotation
alpha = 5*rand(1,1);
beta = 2*rand(1,1);
gamma = 3*rand(1,1);
EditedImage2 = rand(256,1);

[rotatedReal, rotatedImag] = MEX_SOFT_Project(1,8,8,7,alpha, beta, gamma,EditedImage2, zeros(size(EditedImage2)));

fprintf('\nThe result from the new configuration \n')
[new_alpha, new_beta, new_gamma] = MEX_SOFT_Project(0,rotatedReal, rotatedImag,EditedImage2, zeros(size(EditedImage2)), 8,8,7);

disp('These are the actual')
alpha
beta
gamma

disp('These are the estimated')
new_alpha
new_beta
new_gamma


%% Sphere Plot study
N = 200;
deg = 17;
greyColor = [.7 .7 .7];
I = imresize( imadjust(rgb2gray(rot90(imread('world.200412.3x5400x2700.jpg'),2))), [2*N,2*N]);
display_I = double(I);
I = reshape(display_I, [4*N^2,1]);

[x,y,z]= sphere(30);
figure,
hold on
warp(x,y,z,imcomplement(uint8(display_I)));
surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
hold off
axis equal;
grid on
xlabel('x')
ylabel('y')
zlabel('z')

alpha = 5*rand(1,1);
beta = 2*rand(1,1);
gamma = 3*rand(1,1);

[rotatedReal, rotatedImag] = MEX_SOFT_Project(1,N,N,deg,alpha, beta, gamma,I, zeros(size(I)));   %% Rotating syntax

displayRotated = reshape(rotatedReal, [2*N,2*N]);
% displayRotated = uint8(displayRotated);
figure,
hold on
warp(x,y,z,imcomplement(uint8(displayRotated)));
surface(x,y,z,'FaceColor', 'none','EdgeColor',greyColor);
hold off
axis equal;
grid on
xlabel('x')
ylabel('y')
zlabel('z')

fprintf('\nThe result from the new configuration \n')
[new_alpha, new_beta, new_gamma] = MEX_SOFT_Project(0,rotatedReal, rotatedImag,I, zeros(size(I)), N,N,deg);   %% E syntax

disp('These are the actual')
alpha
beta
gamma

disp('These are the estimated')
new_alpha
new_beta
new_gamma
