
%
% create spherical harmonic descriptor
%

function [fvec, deg, Z, new_name] = create_SPHARM_des_LSF(vertices, faces, sph_verts ,maxDeg, filename, OutDirectory)

if isempty(vertices) | isempty(sph_verts)
    if ~isempty(filename)
        load(filename);
    else
        disp('There is no useful information');
        return;
    end
end

[pa, na, ex] = fileparts(filename);
new_name = '';

if ~exist('vertices', 'var') | ~exist('sph_verts', 'var')
    disp('One or more of vertices, and spherical vertices are missing');
    return;
end

%vertnum = size(sph_verts,1);
deg = maxDeg;
%max_d = maxDeg;
% Note that degree 'd' we want to use depends on the vertnum 
% The total number of unknowns is (d+1)*(d+1)
% The total number of equations is vertnum
% We want equ_num >= unk_num
%deg = max(1, floor(sqrt(vertnum)*1/2));
%deg = min(deg, max_d);
% disp(sprintf('Use spharm up to %d degree (vec_len=%d).',deg,(deg+1)^2));

figure, scatter3(vertices(:,1),vertices(:,2),vertices(:,3),'.')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

figure, scatter3(sph_verts(:,1),sph_verts(:,2),sph_verts(:,3),'.')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')


Z = calculate_SPHARM_basis(sph_verts, maxDeg); 

[x,y] = size(Z);
% disp(sprintf('Least square for %d equations and %d unknowns',x,y));

% Least square fitting
% fvec = Z\vertices;   %This does not work as it is expected to work in certain environment
% for i=1:size(vertices,2)
%     fvec(:,i) = Z\vertices(:,i);
% end
 
%//////////////////////////////////////////////////////////////////////////

for i=1:size(vertices,2)
    fvec(:,i) =((Z\vertices(:,i)));
   % fvec(:,i) = (abs(Z\vertices(:,i))).^2;
%      fvec(:,i) =(abs(Z\vertices(:,i)));
end
%Reconstruct
%              lb = 1;
%             ub = (maxDeg+1)^2;
%             verticesR = real(Z(:,lb:ub)*fvec(lb:ub,:));
           % save(objs{i}, 'fvec','vertices','faces', 'sph_verts' );





%  veci=zeros(deg+1,3); %/// change
% %R=zeros(deg,3); %/// change
% for i=1:3
%    
% veci(1,i)=norm((fvec(1,i)));  %1
% veci(2,i)=norm((fvec(2:4,i)));   %3
% veci(3,i)=norm((fvec(5:9,i)));  %5
% veci(4,i)=norm((fvec(10:16,i)));  %7
% veci(5,i)=norm((fvec(17:25,i)));  %9
% veci(6,i)=norm((fvec(26:36,i)));  %11
% veci(7,i)=norm((fvec(37:49,i)));  %13
% veci(8,i)=norm((fvec(50:64,i)));  %15
% veci(9,i)=norm((fvec(65:81,i)));  %17
% veci(10,i)=norm((fvec(82:100,i)));  %19
% veci(11,i)=norm((fvec(101:121,i)));  %21
% veci(12,i)=norm((fvec(122:144,i)));  %23
% veci(13,i)=norm((fvec(144:168,i)));  %25
% veci(14,i)=norm((fvec(169:195,i)));  %27
% veci(15,i)=norm((fvec(196:224,i)));  %29
% veci(16,i)=norm((fvec(225:255,i)));  %31

% vec(1,i)=sqrt(sum(fvec(1,i)));  %1
% vec(2,i)=sqrt(sum(fvec(2:4,i)));   %3
% vec(3,i)=sqrt(sum(fvec(5:9,i)));  %5
% vec(4,i)=sqrt(sum(fvec(10:16,i)));  %7
% vec(5,i)=sqrt(sum(fvec(17:25,i)));  %9
% vec(6,i)=sqrt(sum(fvec(26:36,i)));  %11
% vec(7,i)=sqrt(sum(fvec(37:49,i)));  %13
% vec(8,i)=sqrt(sum(fvec(50:64,i)));  %15
% vec(9,i)=sqrt(sum(fvec(65:81,i)));  %17
% vec(10,i)=sqrt(sum(fvec(82:100,i)));  %19
% vec(11,i)=sqrt(sum(fvec(101:121,i)));  %21
% vec(12,i)=sqrt(sum(fvec(122:144,i)));  %23
% vec(13,i)=sqrt(sum(fvec(144:168,i)));  %25
% vec(14,i)=sqrt(sum(fvec(169:195,i)));  %27
% vec(15,i)=sqrt(sum(fvec(196:224,i)));  %29
% vec(16,i)=sqrt(sum(fvec(225:255,i)));  %31

%end

%Difference
% load ('fvecc.mat');
% for i=1:3
%    
% R(1,i)= veci(1,i)-((vec(1,i)));  %1
% R(2,i)=veci(2,i)-((vec(2,i)));   %3
% R(3,i)=veci(3,i)-((vec(3,i)));  %5
% R(4,i)=veci(4,i)-((vec(4,i)));  %7
% R(5,i)=veci(5,i)-((vec(5,i)));  %9
% R(6,i)=veci(6,i)-((vec(6,i)));  %11
% R(7,i)=veci(7,i)-((vec(7,i)));  %13
% R(8,i)=veci(8,i)-((vec(8,i)));  %15
% R(9,i)=veci(9,i)-((vec(9,i)));  %17
% R(10,i)=veci(10,i)-((vec(10,i)));  %19
% R(11,i)=veci(11,i)-((vec(11,i)));  %21
% R(12,i)=veci(12,i)-((vec(12,i)));  %23
% R(13,i)=veci(13,i)-((vec(13,i)));  %25
% R(14,i)=veci(14,i)-((vec(14,i)));  %27
% R(15,i)=veci(15,i)-((vec(15,i)));  %29
% R(16,i)=veci(16,i)-((vec(16,i)));  %31
% end





%% Save the file
if ~isempty(filename)
    if ~isempty(OutDirectory)
        new_name = sprintf('%s/%s.mat', OutDirectory,  na(1:end)); %filename); %
    else
        new_name = sprintf('%s/%s.mat', pa,na(1:end));    % filename);
    end
%     if exist(new_name, 'file')
%         prompt = {'Enter new filename:'};
%         dlg_title = 'New File Name';
%         num_lines = 1;
%         def = {new_name};
%         answer = inputdlg(prompt,dlg_title,num_lines,def);    
%         new_name = answer{1};
%     end
    save(new_name,'maxDeg', 'vertices', 'faces', 'sph_verts', 'fvec');
end

return;
