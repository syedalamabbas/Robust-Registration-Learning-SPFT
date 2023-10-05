function [dd,thet] = imagecast(data,thetas,nn)
% s=[phi-pi;theta-pi/2; ]';

Coce=1;
data=data';
p=data;

ntheta=360/nn;
nphi=180/nn;
xedge = (0:ntheta:360)*pi/180;
yedge = (0:nphi:180)*pi/180; %linspace(min(data(:,2)),max(data(:,2)),bins);
zedge = linspace(0,2,Coce+1);

xedge=xedge-pi;
yedge=yedge-pi/2;


loc = zeros(size(data));
len1 = length(xedge)-1;
len2 = length(yedge)-1;
len3 = length(zedge)-1;

[~,loc(:,1)] = histc(p(:,1),xedge);
[~,loc(:,2)] = histc(p(:,2),yedge);
[~,loc(:,3)] = histc(p(:,3),zedge);
hasdata = all(loc>0,2);
sz(1:3) = [len1 len2 len3];
loc1=loc(hasdata,:);
thetas=thetas(hasdata);
data=data(hasdata,:);
% cnt = accumarray(loc(hasdata,:),1,sz);
%cnt = accumarray(loc(hasdata,:),1,sz,@mean);

sz=zeros(sz);
cnt1=sz; dd=sz;
loc1(loc1(:,3)>=len3+1,3)=len3;
for i=1:length(thetas)
    sz(loc1(i,1),loc1(i,2),loc1(i,3))=sz(loc1(i,1),loc1(i,2),loc1(i,3))+thetas(i);
    cnt1(loc1(i,1),loc1(i,2),loc1(i,3))=cnt1(loc1(i,1),loc1(i,2),loc1(i,3))+1;
    dd(loc1(i,1),loc1(i,2),loc1(i,3))=max(dd(loc1(i,1),loc1(i,2),loc1(i,3)),data(i,3));
    
end
cnt=cnt1;
sz=sz./cnt;
sz(isnan(sz))=0;
thet=sz;
end


data=data';
p=data;

xedge = linspace(-1,1,nn+1);
yedge = linspace(-1,1,nn+1); %linspace(min(data(:,2)),max(data(:,2)),bins);
zedge = linspace(-1,1,Coce+1);


loc = zeros(size(data));
len1 = length(xedge)-1;
len2 = length(yedge)-1;
len3 = length(zedge)-1;

[~,loc(:,1)] = histc(p(:,1),xedge);
[~,loc(:,2)] = histc(p(:,2),yedge);
[~,loc(:,3)] = histc(p(:,3),zedge);
hasdata = all(loc>0,2);
sz(1:3) = [len1 len2 len3];
%  loc=loc(hasdata,:);
%  datas=data(hasdata);
loc(loc(:,1)>nn,1)=nn;loc(loc(:,2)>nn,2)=nn;loc(loc(:,3)>nn,3)=nn;

cnt = accumarray(loc(hasdata,:),1,sz);
%cnt = accumarray(loc(hasdata,:),1,sz,@mean);
% sz=zeros(sz);
%  cnt1=sz;
%  loc1(loc1(:,3)>=len3+1,3)=len3;
% for i=1:length(datas)
% % sz(loc1(i,1),loc1(i,2),loc1(i,3))=sz(loc1(i,1),loc1(i,2),loc1(i,3))+thetas(i);
% cnt1(loc1(i,1),loc1(i,2),loc1(i,3))=cnt1(loc1(i,1),loc1(i,2),loc1(i,3))+1;
% end
%  cnt=cnt1;
%  cnt=cnt/max(cnt(:));
% sz=sz./cnt;
% sz(isnan(sz))=0;
end