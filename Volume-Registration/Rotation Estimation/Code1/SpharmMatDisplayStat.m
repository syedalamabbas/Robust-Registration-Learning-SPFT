%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayStat(confs, objs, method, windowHandle, cpath)

load(deblank(char(objs)));

[p,fname,e] = fileparts(deblank(char(objs)));

switch lower(method)
    case 'res_t_map'
        switch deblank(char(confs.Overlay))
            case 'p-value'
                cutoff = confs.Threshold_p_value;
                signal = pvalue;
            case 't-map'
                df = length(grInfo)-1;
                cutoff = tinv(1-(confs.Threshold_p_value/2),df);
                signal = tstats;
        end
        patch_overlay(atlas_vertices, faces, signal, cutoff, deblank(char(confs.Overlay)), vtnorm, deblank(char(confs.Colormap)),windowHandle,fname);
    case 'res_pca'
        project_PCA(fvecs, eigenvecs, eigenvals, deblank(char(confs.Mesh)), confs.MaxSPHARMDegree, confs.Level, confs.Sigma, windowHandle,fname, cpath);
end

return;
