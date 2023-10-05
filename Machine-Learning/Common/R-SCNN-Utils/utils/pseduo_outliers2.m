function data=pseduo_outliers2(data,outlayerlevel,range,range2) 




for j=1:length(data(1,1,:))
    
All_points=data(:,:,j)';

n=length(All_points);
nn=fix(outlayerlevel*n); choose=randi(3,1);


for

gnererate    
    



end


All_points=All_points';
All_points=All_points(:,1:(n-length(plane)));
All_points=[All_points,plane'];


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
















