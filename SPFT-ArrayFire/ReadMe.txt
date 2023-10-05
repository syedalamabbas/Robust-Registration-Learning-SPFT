Download and install arrayfire and all of its dependencies from here:
https://arrayfire.com/binaries/

My setup on windows needs these installations:
cmake-3.24.2-windows-x86_64.msi
cuda_11.4.4_472.50_windows.exe
ArrayFire-v3.8.2-CUDA-11.4.exe


To build (the arrayfire based SPFT project) use cmake:

cmake source code directory: SPFT-ArrayFire
cmake source code build directory: SPFT-ArrayFire/build

Configure - Generate - Build 


After it is all done, switch the C++ language feature in general properties
to use C++17