function [pc, score, latent, tsquare] = princomp2(x);
%PRINCOMP Principal Component Analysis (centered and scaled data).
%   [PC, SCORE, LATENT, TSQUARE] = PRINCOMP(X) takes a data matrix X and
%   returns the principal components in PC, the so-called Z-scores in SCORES,
%   the eigenvalues of the covariance matrix of X in LATENT, and Hotelling's
%   T-squared statistic for each data point in TSQUARE.

%   Reference: J. Edward Jackson, A User's Guide to Principal Components
%   John Wiley & Sons, Inc. 1991 pp. 1-25.

%   B. Jones 3-17-94
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.9 $  $Date: 2002/01/17 21:31:45 $

[m,n] = size(x);
r = min(m-1,n);     % max possible rank of x
avg = mean(x);
centerx = (x - avg(ones(m,1),:));

[U,latent,pc] = svd(centerx./sqrt(m-1),0);
score = centerx*pc;

if nargout < 3, return; end
latent = diag(latent).^2;
if (r<n)
   latent = [latent(1:r); zeros(n-r,1)];
   score(:,r+1:end) = 0;
end

if nargout < 4, return; end
tmp = sqrt(diag(1./latent(1:r)))*score(:,1:r)';
tsquare = sum(tmp.*tmp)';
