%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% fix  vertices
%

function [vec1, vec2] = fixed_vec(vec1,vec2)


len =size(vec1,1);

if size(vec2,1)>=len
    vec2 = vec2(1:len,:);
else
   len =size(vec2,1);
   vec1 = vec1(1:len,:);
end
   
return;