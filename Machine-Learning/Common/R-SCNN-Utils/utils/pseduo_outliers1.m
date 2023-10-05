function data=pseduo_outliers1(data,outlayerlevel,range) 


for j=1:length(data(1,1,:))
    
All_points=data(:,:,j)';

n=length(All_points);
nn=fix(outlayerlevel*n); choose=randi(3,1);

%p=data;
if choose==1
% genrerate plane in z
D=min(All_points(:,3));
D=D+range*D/abs(D);
%A=0;B=0;C=1;
x=min(All_points(:,1))+(max(All_points(:,1))-min(All_points(:,1))).*rand(nn,1);
y=min(All_points(:,2))+(max(All_points(:,2))-min(All_points(:,2))).*rand(nn,1);
%z = -1/C*(A*x + B*y + D);
z=ones(nn,1)*D;
plane=[x,y,z];

elseif choose==2
% genrerate plane in x
D=min(All_points(:,1));
D=D+range*D/abs(D);
z=min(All_points(:,3))+(max(All_points(:,3))-min(All_points(:,3))).*rand(nn,1);
y=min(All_points(:,2))+(max(All_points(:,2))-min(All_points(:,2))).*rand(nn,1);
x  = ones(nn,1)*D;
plane=[x,y,z];

else
   
% genrerate plane in y
D=min(All_points(:,2));
D=D+range*D/abs(D);

z=min(All_points(:,3))+(max(All_points(:,3))-min(All_points(:,3))).*rand(nn,1);
x=min(All_points(:,1))+(max(All_points(:,1))-min(All_points(:,1))).*rand(nn,1);
y  = ones(nn,1)*D;
plane=[x,y,z];
  
end
All_points=All_points';
All_points=All_points(:,1:(n-length(plane)));
All_points=[All_points,plane'];

% scatter3(x,y,z);
% hold on
% plot3(p(:,1),p(:,2),p(:,3),'g.')
% xlabel('x') 
% ylabel('y') 
% zlabel('z')     
%  hold off  


% i=randperm(length(All_points));
% All_points=All_points(:,i);
% All_points=All_points(:,1:n);

%% save data to same file
data(:,:,j)=All_points;

% scatter3(All_points(1,:),All_points(2,:),All_points(3,:))
%         xlabel('My x label')
%         ylabel('y')
%         zlabel('zz')
end