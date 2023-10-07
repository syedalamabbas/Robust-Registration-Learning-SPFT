function  SpecialIsoSurfaceDisplay( volume_3D, surfaceVal )
%SPECIALISOSURFACEDISPLAY Summary of this function goes here
figure,
axis tight
hold on
hpatch = patch(isosurface(volume_3D,surfaceVal));
isonormals(volume_3D,hpatch)
hpatch.FaceColor = 'cyan';
hpatch.EdgeColor = 'none';
daspect([1,1,1.7])
view([-57,68])
camlight left;
lighting gouraud
hold off
grid on
xlabel('x')
ylabel('y')
zlabel('z')

end

