function [ normals ,D] = findPointNormals(points, numNeighbours, viewPoint, dirLargest)
%FINDPOINTNORMALS Estimates the normals of a sparse set of n 3d points by
% using a set of the closest neighbours to approximate a plane.
%
%   Required Inputs:
%   points- nx3 set of 3d points (x,y,z)
%
%   Optional Inputs: (will give default values on empty array [])
%   numNeighbours- number of neighbouring points to use in plane fitting
%       (default 9)
%   viewPoint- location all normals will point towards (default [0,0,0])
%   dirLargest- use only the largest component of the normal in determining
%       its direction wrt the viewPoint (generally provides a more stable
%       estimation of planes near the viewPoint, default true)
%
%   Outputs:
%   normals- nx3 set of normals (nx,ny,nz)
%   curvature- nx1 set giving the curvature
%
%   References-
%   The implementation closely follows the method given at
%   http://pointclouds.org/documentation/tutorials/normal_estimation.php
%   This code was used in generating the results for the journal paper
%   Multi-modal sensor calibration using a gradient orientation measure 
%   http://www.zjtaylor.com/welcome/download_pdf?pdf=JFR2013.pdf
%
%   This code was written by Zachary Taylor
%   zacharyjeremytaylor@gmail.com
%   http://www.zjtaylor.com

%% check inputs
validateattributes(points, {'numeric'},{'ncols',3});

if(nargin < 2)
    numNeighbours = [];
end
if(isempty(numNeighbours))
    numNeighbours = 9;
else
    validateattributes(numNeighbours, {'numeric'},{'scalar','positive'});
    if(numNeighbours > 100)
        warning(['%i neighbouring points will be used in plane'...
            ' estimation, expect long run times, large ram usage and'...
            ' poor results near edges'],numNeighbours);
    end
end

if(nargin < 3)
    viewPoint = [];
end
if(isempty(viewPoint))
    viewPoint = [0,0,0];
else
    validateattributes(viewPoint, {'numeric'},{'size',[1,3]});
end

if(nargin < 4)
    dirLargest = [];
end
if(isempty(dirLargest))
    dirLargest = true;
else
    validateattributes(dirLargest, {'logical'},{'scalar'});
end

%% setup

%ensure inputs of correct type
points = double(points);
viewPoint = double(viewPoint);

%create kdtree
kdtreeobj = KDTreeSearcher(points,'distance','euclidean');

%get nearest neighbours
n = knnsearch(kdtreeobj,points,'k',(numNeighbours+1));

%remove self
n = n(:,2:end);

%find difference in position from neighbouring points
p = repmat(points(:,1:3),numNeighbours,1) - points(n(:),1:3);
p = reshape(p, size(points,1),numNeighbours,3);

%calculate values for covariance matrix
C = zeros(size(points,1),6);
C(:,1) = sum(p(:,:,1).*p(:,:,1),2);
C(:,2) = sum(p(:,:,1).*p(:,:,2),2);
C(:,3) = sum(p(:,:,1).*p(:,:,3),2);
C(:,4) = sum(p(:,:,2).*p(:,:,2),2);
C(:,5) = sum(p(:,:,2).*p(:,:,3),2);
C(:,6) = sum(p(:,:,3).*p(:,:,3),2);
C = C ./ numNeighbours;

%% normals and curvature calculation

normals = zeros(size(points));
D = zeros(size(points,1),1);
% curvature = zeros(size(points,1),1);
for i = 1:(size(points,1))
    
    %form covariance matrix
    Cmat = [C(i,1) C(i,2) C(i,3);...
        C(i,2) C(i,4) C(i,5);...
        C(i,3) C(i,5) C(i,6)];  
    
    %get eigen values and vectors
    [v,d] = eig(Cmat);
    d = diag(d);
    [lambda,k] = min(d);
    
    %store normals
    normals(i,:) = v(:,k)';
%     D(i)=(normals(i,1)*(points(i,1))+normals(i,2)*(points(i,2))+normals(i,3)*(points(i,3)))/ sqrt(normals(i,1)^2+normals(i,2)^2+normals(i,3)^2);
%     theta1(i,:)=thetas(normals(i,:))';
% %     yaw = theta(i,1);
% %     pitch = theta(i,2);
% %     roll = theta(i,3);
%     
% %         X = sin(yaw);
% %         Y = -(sin(pitch)*cos(yaw));
% %         Z = -(cos(pitch)*cos(yaw));
% %      r(i,:)=[X Y Z];
%      v=[points(i,:)' points(n(i,:),:)'];
%      [inl, plane(i,:)]=ransaymanold(v,.01,30);
%     phi= theta1(i,3);
%     theta=theta1(i,1);
%  beta=cos(deg2rad(phi))*points(i,3)+ sin(deg2rad(phi))*cos(deg2rad(theta))*points(i,1)+sin(deg2rad(phi))*sin(deg2rad(theta))*points(i,2);  
% D(i)=beta;
    %     [yaw, pitch, roll] = rod2angle(normals(i,:))
    %store curvature
%     curvature(i) = lambda / sum(d);
end

%% flipping normals

%ensure normals point towards viewPoint
% points = points - repmat(viewPoint,size(points,1),1);
% if(dirLargest)
%     [~,idx] = max(abs(normals),[],2);
%     idx = (1:size(normals,1))' + (idx-1)*size(normals,1);
%     dir = normals(idx).*points(idx) > 0;
% else
%     dir = sum(normals.*points,2) > 0;
% end
% 
% normals(dir,:) = -normals(dir,:);
% % D=(normals(:,1).*(points(:,1))+normals(:,2).*(points(:,2))+normals(:,3).*(points(:,3)));%./ sqrt(normals(:,1).^2+normals(:,2).^2+normals(:,3).^2);
% 
% % rod=angle2rod(theta(1,:))
% 
% % r = angle2rod(yaw,pitch,roll);
% end
% viewPoint=[0 0 10];
% points = points - repmat(viewPoint,size(points,1),1);
% if(true)
%     [~,idx] = max(abs(normals),[],2);
%     idx = (1:size(norma1,1))' + (idx-1)*size(norma1,1);
%     dir = norma1(idx).*points(idx) > 0;
% else
%     dir = sum(norma1.*points,2) > 0;
% end
% 
% norma1(dir,:) = -norma1(dir,:);

