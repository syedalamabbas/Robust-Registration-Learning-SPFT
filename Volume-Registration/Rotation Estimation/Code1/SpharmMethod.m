

function [rmsd, theta, R] = SpharmMethod(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end
outNames = {};

%class(confs) - struct
%confs.vars -list of variable
%class(objs) - cell

% h = waitbar(0,'Please wait...');

%% one or many fils
%%
 for i = 1: numSbj
%  for i = 1: 1
    file =  objs{i};
%     file =  objs;
    [path, name, ext] = fileparts(file);
   switch method
        case 'Eigenvector'
            % Call Eigenvector method
            if isempty(confs.Template)
                disp('Eigenvector needs a template object');
                return;
            end
               [rmsd(1,i),theta, R] = eigenvector(file, confs);
%              [rmsd(1,i),theta, R]=  eigenvector2(file, confs);
                 %break;
        case 'FOE'
            % Call FOE method
             [rmsd(1,i),theta,R] = FOE(file, confs);
%             [rmsd(1,i),theta,R] = align_FOE2(file, confs);     
            % break;
                 
       case 'ICP'
            % Call Eigenvector method
            if isempty(confs.Template)
                disp('ICP_Eigenvector needs a template object');
                return;
            end
%          [rmsd(1,i),theta,R] = ICP_eigenvector(file, confs); % using fvec  
 %         [rmsd(1,i),theta,R]= ICP1(file, confs);   % using vertices iterative
           [rmsd(1,i),theta,R]  = ICP2(file, confs); % using{'bruteForce','normalShooting','Delaunay','GLtree','kDtree'};


                % break;
       case 'PCA'
            [rmsd(1,i),theta,R] = pca(file, confs); 
           % break;
            
       otherwise
              disp('No method defined.');
              rmsd=0;
    end  
    
%       diary('off');
%   waitbar(i/numSbj)
end
 %close(h);

return