function hFigOriginal = Visualize3SlicePlanesOfAVolume( simpleVolume )
%VISUALIZE3SLICEPLANESOFAVOLUME displays 3 slices at the origin from a
%volume

sizeIn = size(simpleVolume);
hFigOriginal = figure;
hAxOriginal  = axes;
slice(double(simpleVolume),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
xlabel('x')
ylabel('y') 
zlabel('z')
grid on,
shading interp,
colormap hot
rotate3d on
axis equal

end

