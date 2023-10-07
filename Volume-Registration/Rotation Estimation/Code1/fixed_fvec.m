%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% fix  fvec
%

function [fvec, max_d] = fixed_fvec(fvec, max_d, scale)

% fvec(1,:) = 0;
len = (max_d+1)^2;
% len1 = (11)^2;
if size(fvec,1)>len
    fvec = fvec(1:len,:);
% else
%     disp('lenght of fvec is less than sphram band');
end
max_d = min(max_d,sqrt(size(fvec,1))-1);


    
return;