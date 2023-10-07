function [ translation_vector ] = ComputeTranslationWithVolumes( volume_original, volume_translated )
%COMPUTETRANSLATIONWITHVOLUMES function computes the translation as
%described in the paper IEEE PAMI by Syed Alam Abbas
%=====================================================================%
% Written on 6/22/2016
% Updated by adding function ComputeFreq1DData on 7/8/2016
% Revised on March 12th, 2022 by Syed Alam Abbas.
%=====================================================================%

% [N, ~, ~] = size(volume_original);
% if(mod(N,2))  % Precautionary measure for odd case
%     volume_original = padarray(volume_original, [1 1 1],'pre');
%     volume_translated = padarray(volume_translated, [1 1 1],'pre');
% end
% [N, ~, ~] = size(volume_original);  % Definitely is assumed to be N = even

%% Taking FFT
FFT_original = fftn(volume_original);             % 3D FFT
FFT_translated = fftn(volume_translated);         % 3D FFT
PhaseCorrelation = (FFT_translated.* conj(FFT_original))./ abs(FFT_translated .*conj(FFT_original));        % According to the formula given in the paper

InversePhaseCorrelation = ifftshift(ifftn(PhaseCorrelation) );
[m,ix] = max(InversePhaseCorrelation(:));
[tx,ty,tz]= ind2sub(size(InversePhaseCorrelation),ix);
translation_vector  = [tx, ty, tz];
disp(['Computed translation vector=[', num2str(tx), ',', num2str(ty),',', num2str(tz), ']' ]);

% %% Computing the translations now from the Phase Correlations using our refined approach
% % Line_index = N/2 + 5;                    % This is an arbitrary index , gotta change this to cover the entire matrix of Phase Correlation
% 
% %% X axis
% LineData_Row = real(PhaseCorrelation( :, ty, tz));                 % Just the cosine part that has the phase information encoded as translation
% est_tx  = ComputeFreq1DData( LineData_Row, N );
% 
% %% Y axis
% LineData_Column = real(PhaseCorrelation( tx,:,tz));
% est_ty  =  ComputeFreq1DData( LineData_Column, N );
% 
% %% Z axis
% LineData_Depth = real(PhaseCorrelation( tx,ty,:));
% est_tz  = ComputeFreq1DData( LineData_Depth, N );
% 
% 
% %% Final answer - Optimal translation
% disp(['Computed translation vector refined =[', num2str(est_tx), ',', num2str(est_ty),',', num2str(est_tz), ']' ]);
% 
%     function FreqEstimate = ComputeFreq1DData(LineData , N)
%         resolution_factor = 128;                                   % Can change it to arbitrary precision
%         [S,w]  = pmusic(LineData(1:N/2),2, N*resolution_factor);   % Using the MATLAB pmusic method to compute the frequency of the 1-D signal
%         w = w/pi;
%         maxIndex = (find(S == max(S)));
%         f_estimate1 = w(maxIndex);
%         FreqEstimate  = f_estimate1*N;
%     end
end

