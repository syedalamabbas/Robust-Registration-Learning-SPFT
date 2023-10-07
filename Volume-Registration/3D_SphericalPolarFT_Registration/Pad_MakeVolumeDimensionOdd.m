function [ final_volume1, final_volume2 , N  ] = Pad_MakeVolumeDimensionOdd( volume1, volume2 )
%PAD_MAKEVOLUMEDIMENSIONODD 
% This function does central padding of the volumes 1 and 2 such that the
% dimensions of both of them become equal 
% It also ensures that the final volumes 1 and 2 are of N^3 size each where
% N = odd
%====================================================================
% Written on June 13th, 2016 by Syed Alam Abbas.
% Revised on March 12th, 2022 by Syed Alam Abbas.
%====================================================================

[ v1_x, v1_y, v1_z ] = size(volume1); 
[ v2_x, v2_y, v2_z ] = size(volume2);

%% Make all the dimensions odd - all may still be different
final_volume1 = volume1;        % Initilize
final_volume2 = volume2;

if(~mod(v1_x,2))  % Not an odd number
    final_volume1 = padarray(final_volume1, [1 0 0],'pre');
end
if(~mod(v2_x,2))  % Not an odd number
    final_volume2 = padarray(final_volume2, [1 0 0],'pre');
end

if(~mod(v1_y,2))  % Not an odd number
    final_volume1 = padarray(final_volume1, [0 1 0],'pre');
end
if(~mod(v2_y,2))  % Not an odd number
    final_volume2 = padarray(final_volume2, [0 1 0],'pre');
end

if(~mod(v1_z,2))  % Not an odd number
    final_volume1 = padarray(final_volume1, [0 0 1],'pre');
end
if(~mod(v2_z,2))  % Not an odd number
    final_volume2 = padarray(final_volume2, [0 0 1],'pre');
end


%% Now we can equalize odd dimensions volumes using max values
[ v1_x, v1_y, v1_z ] = size(final_volume1); 
[ v2_x, v2_y, v2_z ] = size(final_volume2);

v2_x_max = max(v1_x, v2_x);
v2_y_max = max(v1_y, v2_y);
v2_z_max = max(v1_z, v2_z);

%% X - axis equalization
if(v2_x_max > v1_x)
    final_volume1 = padarray(final_volume1, [(v2_x_max-v1_x)/2 0 0],'both');
else
    final_volume2 = padarray(final_volume2, [(v2_x_max-v2_x)/2 0 0],'both');
end

%% Y - axis equalization
if(v2_y_max > v1_y)
    final_volume1 = padarray(final_volume1, [ 0 (v2_y_max-v1_y)/2 0],'both');
else
    final_volume2 = padarray(final_volume2, [ 0 (v2_y_max-v2_y)/2 0],'both');
end

%% Z - axis equalization
if(v2_z_max > v1_z)
    final_volume1 = padarray(final_volume1, [0 0 (v2_z_max-v1_z)/2 ],'both');
else
    final_volume2 = padarray(final_volume2, [ 0 0 (v2_z_max-v2_z)/2 ],'both');
end

%% Final Equalization for N^3 volumes
[ v1_x, v1_y, v1_z ] = size(final_volume1); 
[ v2_x, v2_y, v2_z ] = size(final_volume2);

N = max( [v1_x, v1_y, v1_z, v2_x, v2_y, v2_z] );

%% X - axis equalization
if(N > v1_x)
    final_volume1 = padarray(final_volume1, [(N-v1_x)/2 0 0],'both');
end
if(N > v2_x)
    final_volume2 = padarray(final_volume2, [(N-v2_x)/2 0 0],'both');
end

%% Y - axis equalization
if(N > v1_y)
    final_volume1 = padarray(final_volume1, [ 0 (N-v1_y)/2 0],'both');
end
if(N > v2_y)
    final_volume2 = padarray(final_volume2, [ 0 (N-v2_y)/2 0],'both');
end

%% Z - axis equalization
if(N > v1_z)
    final_volume1 = padarray(final_volume1, [0 0 (N-v1_z)/2 ],'both');
end
if(N > v2_z)
    final_volume2 = padarray(final_volume2, [ 0 0 (N-v2_z)/2 ],'both');
end

end

