function [data]=mergedata(x1,x2)
%Merging in rows, columns should be same
%
%
if(x1==0)
    data=x2;
else
[row1,col1]=size(x1);
[row2,col2]=size(x2);
if(col1~=col2)
    error('Two matrix columns are not same', ...
      'There exists an invalid input.');
end
data=zeros(row1+row2,col1);
data(1:row1,1:col1)=x1(:,:);
data(row1+1:(row1+row2),1:col1)=x2(:,:);
end


