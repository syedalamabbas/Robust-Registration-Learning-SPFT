
close all
clc
clear


%% Random vector color plot

data = rand(1,100); %// random dummy data - 100 element vector
im = imagesc(data);
colorbar;
set(gca, 'YTick', []);


% Make some 4D data in X, Y, height, colour:
% Each point now has an X, Y, Z (height), and C (colour) value
[X,Y,Z] = peaks(51);    % 51-by-51 Z as a function of X, Y
C = membrane(1,25);     % Colour
% Show the data as a surface plot
figure, surf(X,Y,Z,C), colorbar
% Show the data as a scatter3 plot
X_flat = X(:);
figure, scatter3(X(:),Y(:),Z(:),6,C(:),'filled'), colorbar


