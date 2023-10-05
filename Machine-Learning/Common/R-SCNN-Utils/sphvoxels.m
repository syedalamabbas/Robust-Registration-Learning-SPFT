function [cnt,sz] = sphvoxels(data,thetas,nn,Coce)
% s=[phi-pi;theta-pi/2; ]';

data=data';
p=data;

ntheta=360/nn;
nphi=180/nn;

ntheta1=360/(nn-1);
nphi1=180/(nn-1);
xedge = (ntheta/2:ntheta1:(360+ntheta/2))*pi/180; 
yedge = (nphi/2:nphi1:(180+nphi/2))*pi/180; 

zedge = linspace(0,1.5,Coce+1);
zedge=zedge+zedge(2)/2;
zedge=zedge(1:end-1);

xedge=xedge-pi;
yedge=yedge-pi/2;

loc = zeros(size(data));
len1 = length(xedge);
len2 = length(yedge);
len3 = length(zedge);

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
 cnt1=sz;
 loc1(loc1(:,3)>=len3+1,3)=len3;
for i=1:length(thetas)
    if(loc1(i)<=nn)
 sz(loc1(i,1),loc1(i,2),loc1(i,3))=sz(loc1(i,1),loc1(i,2),loc1(i,3))+thetas(i);
cnt1(loc1(i,1),loc1(i,2),loc1(i,3))=cnt1(loc1(i,1),loc1(i,2),loc1(i,3))+1;
    end
end
 cnt=cnt1;
  sz=sz./cnt;
  sz(isnan(sz))=0;
end