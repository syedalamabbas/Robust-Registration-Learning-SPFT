function [ x_Solution ] = SymbolicComputationFactorizeRotation( knownT )
%% Given Transformation
% knownT = 3 x 3 transformation
% knownT = [ 0.9704   -0.0143   -0.2410 ;
%     0.0228    0.9992    0.0324 ;
%     0.2404   -0.0369    0.9700 ];

syms alpha beta gamma

localTestFun = @(x) GetFullRotationMatrixZYZ( x(1), x(2), x(3)) - knownT ;

x0 = [0.6,0.7,0.8];
x_Solution = fsolve(localTestFun,x0);

disp(x_Solution)
end
