# superscs_win
Python module using superscs on WindowsX64

This python module is to call superscs optimizer on WindowsX64, with help of nimporter and nimpy.
Python calls a compiled extension module written in Nim, which in turn loads shared library of superscs.

Somehow superscs is slower than scs, maybe because I'm using Ryzen CPU, and OpenBLAS is not as 
efficient as BLAS-MKL. I don't know why.
