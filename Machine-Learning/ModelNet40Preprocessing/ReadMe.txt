Downloaded from 
https://www.kaggle.com/datasets/balraj98/modelnet40-princeton-3d-object-dataset

Problems with this data set: 

Some OFF files have seperate line definition for vertices and faces,
e.g. airplane
OFF
90714 104773 0

But some OFF files have same line definition,
e.g. bathtub
OFF3514 3546 0  

For this the read_off file has been modified to match this.

Inconsistencies such as: 
tv_stand is used in the path of metadata file but the folder name is tv
flower_pot is flower, glass_box is glass etc.

have to be resolved.

Then the sequence of operations occurs as follows:
1. archive has off files they are converted into logical volumes.
2. extracted_volumes have 65 x 65 x 65 size, and are transformed into spherical grid.
3. extracted_spherical_grids have 64 x 32 x 65 size, and it is magnitude of the Spherical Polar Fourier Transform.
4. extracted_multispectral_sp_images have 64 x 64 x 11 size, they can now be used for training and testing, etc.