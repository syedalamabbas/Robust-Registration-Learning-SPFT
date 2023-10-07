function [Ept]=extreme_points(d)
y=d(find(d~=0));
[row,col]=size(y);
x=1:col;
%plot(x,y,'-');
%grid on;
%hold on;
t=fpeak(x,y,10,[1,col,30,inf]);
%plot(t(:,1),t(:,2),'o');
%title('\fontsize{24}Perfect!');
[row,col]=size(t);
 % Ept(1)=row;
for ii=1:row
  Ept(ii)=t(ii,2);
 % Ept(row+ii)=t(ii,2);
end

%hold off;