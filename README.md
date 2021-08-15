# superscs_win
This is a Python3 module to use [Superscs](https://github.com/kul-optec/scs) (a conic optimizer written in C) on WindowsX64

The original Superscs package has a python setup.py file, but when I try to install it on WinX64, there are some werid errors, which are likely caused by Microsoft compilers (default C compiler used by setuptools).  

Anyway, so I tried a different approach: compile superscs as a shared library, use nim to create a
superscs api wrapper, then use nimporter/nimpy to load superscs nim functions as python module.

This python module is to call superscs optimizer on WindowsX64, with help of nimporter and nimpy.
Python calls a compiled extension module written in Nim, which in turn loads shared library of superscs.

Somehow superscs is slower than scs, maybe because I'm using Ryzen CPU, and OpenBLAS is not as 
efficient as BLAS-MKL. I don't know why. Before knowing why superscs is slower on my computer, this package is just a practice project to write a python extension using Nim.
