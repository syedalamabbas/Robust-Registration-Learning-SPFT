% ============================================
%
% Goal: match P to X 
%


function [P,M] = match_ICP_eigen(P,X)

M = eye(4);

while(1)   
   for k = 1 : length(P)
      d = (X(:,1)-P(k,1)).^2 + (X(:,2)-P(k,2)).^2 + (X(:,3)-P(k,3)).^2;
      [val,ind] = min(d);
      c(k) = ind;
   end
   Xc = X(c,:);
   
   
   [R,qT] = decomp_eigen(P,Xc);

   Pnew  = (R*P')';
   diff  = mean( sum(((Pnew - P).^2)') );


   P = Pnew;
   
   P(:,1) = P(:,1) + qT(1);
   P(:,2) = P(:,2) + qT(2);
   P(:,3) = P(:,3) + qT(3);
   
   Mr = eye(4); Mr(1:3,1:3) = R;
   Mt = eye(4); Mt(1:3,4) = qT(1:3);
   M = Mt*Mr*M;
   if( diff < 0.00000000001 )
    %  fprintf( 'difference=%f\n', diff );
      break;
   end
   
end

