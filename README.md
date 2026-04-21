# WaveletMinicourse

A minicourse to help CPU and FPGA programmers understand, use and implement DWT for JPEG XS codecs.

* In the executable 'haar-int2int', the Haar wavelet transform operates over a ring of integers: 
see comments in lines 29-30 and 61-62 of 'haar-int2int.cpp' for explanation.

* A pair of direct and inverse Le Gall 5/3 transforms as implemented in the 'legall53' folder (and 
similar ones throughout this project) guarantees a lossless coding over a ring of integers. The 
trick used there is summarized in 
[THE LIFTING SCHEME: A CONSTRUCTION OF SECOND GENERATION WAVELETS by Wim Sweldens](https://cm-bell-labs.github.io/who/wim/papers/lift2.pdf) 
on page 35:

`14.11. Integer to integer wavelet transforms. In [16] lifting is used to build reversible
wavelets which map integers to integers for applications to lossless image coding. The idea is
to introduce a non-linear round-off in each lifting step. This way the result is guaranteed to 
be integer while the lifting assure that the transform is invertible. The exact same idea works 
in the second generation setting and using this in second generation compression applications 
is another line of future research.`

