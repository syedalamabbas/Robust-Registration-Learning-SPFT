

%
% match fvec to atlas by rotating the parameterization of fvec
% - object does not transform in the object space
% - alpha, beta, gamma: samples in the rotation space
% - gran = # of alphas to process together
% - N: keep top N results
%

function [best_fvec, best_angl] = match_param(fvec,atlas,alpha,beta,gamma,gran,N,max_d)

% use which degree for SHREC
deg = 1;
deg = deg+1;

best_fvec = fvec; best_fvec = best_fvec(:,:,ones(1,N));

rsmd = SPHARM_rmsd(fvec, atlas);

% rsmdvec = SPHARM_rmsdvec(fvec, atlas);
% rsmd = rsmdvec(deg);

best_rmsd = ones(1,N)*rsmd;
best_angl = zeros(N,3);

% find best correspondence by rotating parameterization
n = length(alpha); m = length(gamma); nm = n*m;
alpha = alpha(:,ones(1,m)); beta = beta(:,ones(1,m)); gamma = gamma(ones(1,n),:);
xx = ceil(nm/gran); bd = [0 (1:xx)*gran]; bd(end) = nm;    
for k = 1:xx
    idx = (bd(k)+1):bd(k+1); g = length(idx);
    [rvec] = rotate_param_m03(alpha(idx)', beta(idx)', gamma(idx)', fvec, max_d); % m03 < m02 < m01 < s03 < s01 < s00

    rmsd = SPHARM_rmsd_batch(rvec, atlas);
    
%     rmsdvec = SPHARM_rmsdvec_batch(rvec, atlas);
%     rmsd = rmsdvec(deg,:);
    
    best_rmsd = [best_rmsd rmsd];
    best_fvec(:,:,end+(1:g)) = rvec;
    best_angl(end+(1:g),:) = [alpha(idx)', beta(idx)', gamma(idx)'];
    [best_rmsd,ind] = sort(best_rmsd); best_rmsd = best_rmsd(1:N); ind = ind(1:N);
    best_fvec = best_fvec(:,:,ind);
    best_angl = best_angl(ind,:);
  %  disp([sprintf('(%d/%d) RMSD = %f, BEST =',bd(k)+j,nm,rmsd(end)), sprintf(' %f',best_rmsd)]);
end
    
return;