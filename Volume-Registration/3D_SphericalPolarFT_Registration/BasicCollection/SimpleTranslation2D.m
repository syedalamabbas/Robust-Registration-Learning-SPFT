clear all;
close all;
clc;

%% Testing just 2D here
Nx =  256;          % Size of the image
t_x = 56.43;        % x location
t_y = 49.43;        % y location

%% Tale of Two images
phase = 0;
fixedI  = imresize( imread('cameraman.tif') , [Nx,Nx]);
ref = imref2d(size(fixedI));
fixedI  = im2double(fixedI  );
[nr,nc]=size(fixedI);
Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
[Nc,Nr] = meshgrid(Nc,Nr);
movingI = ifft2(fft2(fixedI).*exp(1i*2*pi*(t_x*Nr/nr+t_y*Nc/nc))).*exp(-1i*phase);

snr = 60;
movingI = awgn(movingI,snr,'measured') ;  % infecting the signal
figure(1);
subplot(1,2,1);
imshow(abs(fixedI),ref);
title('Reference image, I_1(x,y)')
subplot(1,2,2);
imshow(abs(movingI),ref);
str = strcat(' I_2(x,y): t_x=',num2str(t_x),',t_y =',  num2str(t_y)  );
title(str)

%% Using PAMI External Code by Dr. Vas
addpath(genpath('..\..\2D_PolarFT_Registration\2010PAMI_RobustFFTBasedRegistration'));
I1 = fixedI;
I2 = movingI;
[scale, th, dis]=PAMIvas(I1,I2);

t_y_estimated = dis(1)
t_x_estimated = dis(2)


%% Using my code now 
%% Taking FFT
FFT_fixedI = (fft2(fixedI));
FFT_movingI = (fft2(movingI));
PhaseCorrelation = (FFT_movingI.* conj(FFT_fixedI))./ abs(FFT_movingI .*conj(FFT_fixedI));        % According to the formula given in the paper
figure, imagesc(log(abs(PhaseCorrelation)))

%% This is the expected correlation
nr = Nx;
nc = Nx;
% [nr,nc] = [512, 512];
Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
[u,v] = meshgrid(Nc,Nr); 

% ExpectedCorrelation =  exp(1i*2*pi*(t_x*u/Nx+t_y*v/Nx));    % This is from exact direct solution
ExpectedCorrelation  = PhaseCorrelation;           % At very low noise this technique works

figure, imagesc( (abs(ExpectedCorrelation) ))
colormap(hot) 

figure, imshow( (angle(ExpectedCorrelation) ))
xlabel('u')
ylabel('v')
colormap(hot)
 
Line_index = 144;          % this is an arbitrary index

%% Column operations for computing delta y or t_y
LineData_Column = real(ExpectedCorrelation( Line_index,:));
figure, plot(LineData_Column, 'LineWidth',1.5)
xlabel('v')
ylabel('Intensity')
grid on
 
% n = 0:128;
% y = sin(.2926*pi*n)+.01 * rand(1, length(n));
% figure, pmusic(y,2, 4096)             % Changing NFFT values directly  increase the resolution of pmusic method
resolution_factor = 128;                 % Can change it to arbitrary precision
 
figure,
[S,w]  = pmusic(LineData_Column(1:Nx/2),2, Nx*resolution_factor);
hold on
w = w/pi;
maxIndex = (find(S == max(S))); 
f_estimate1 = w(maxIndex);                 %% 0.2927 is the peak value for t_y = 37.46, check it 
plot (w ,20*log10(abs(S)),w(maxIndex), 20*log10(abs(S (maxIndex))), '*',w(maxIndex), 20*log10(abs(S (maxIndex))), 'o')     % Normalized Frequency vs Pseudo Spectrum 
est_ty  = f_estimate1*Nx/2 
str1 = strcat('\leftarrow t_y = ',num2str(est_ty) );
text(w(maxIndex), 20*log10(abs(S (maxIndex))),str1)

grid on 
% title('Computing the translation t_y using 1D pseudo-spectrum of Cross Correlation')
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Power (dB)')
hold off
 
%% Row operation for computing delta y or t_y
LineData_Row = real(ExpectedCorrelation( :,Line_index));
figure, plot(LineData_Row, 'LineWidth',1.5)
xlabel('u') 
ylabel('Intensity')
grid on  

figure, 
[S,w]  =  pmusic(LineData_Row(1:Nx/2),2,Nx*resolution_factor);
hold on
w = w/pi;
maxIndex = (find(S == max(S)));
f_estimate2 = w(maxIndex);                 %% .1831 is the peak value for t_x = 23.43, you can check it 
plot (w ,20*log10(abs(S)),w(maxIndex), 20*log10(abs(S (maxIndex))), '*',w(maxIndex), 20*log10(abs(S (maxIndex))), 'o')     % Normalized Frequency vs Pseudo Spectrum 
est_tx  = f_estimate2*Nx/2
str1 = strcat('\leftarrow t_x = ',num2str(est_tx) );
text(w(maxIndex), 20*log10(abs(S (maxIndex))),str1)
grid on 
% title('Computing the translation t_x using 1D pseudo-spectrum of Cross Correlation')
xlabel('Normalized Frequency (\times \pi rad/sample)')
ylabel('Power (dB)')
hold off

%% Printing files

% print -dpdf GivenImagesTranslation
% print -dpdf RowData
% print -dpdf RowPeak
% print -dpdf ColumnData
% print -dpdf ColumnPeak
