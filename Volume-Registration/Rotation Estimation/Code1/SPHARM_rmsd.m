

%
% distance between two SPHARM surfaces
%

function rmsd = SPHARM_rmsd(fvec1, fvec2)

dist = fvec1-fvec2;
%rmsd = norm(dist(:));
% rmsd = norm(dist(:))/sqrt((max_d+1)^2);
rmsd = norm(dist(:))/sqrt(4*pi);
return;

