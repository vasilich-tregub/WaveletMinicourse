#include <stdio.h>
#include <vector>
#include "LGT.h"
#include "CDF97.h"

int main() {
    const int n = 64;
    std::vector<double> im53(n);
    std::vector<double> im97(n);
    int i;

    // Makes a polinomial signal
    for (i = 0; i < n; i++) im97[i] = im53[i] = 5 + i + 0.4 * i * i -0.02 * i * i * i;

    // Prints original signal im
    printf("Original signal:\n");
    for (i = 0; i < n; i++) printf("im[%d]=%f\n", i, im53[i]);
    printf("\n");

    // Do the forward transforms
    fwt53(im53, 0);
    fwt97(im97, 0);

    // Prints the wavelet coefficients
    printf("level 0 Wavelets coefficients:\n");
    for (i = 0; i < n; i++) printf("wc[%d] = %f ; %f\n", i, im53[i], im97[i]);
    printf("\n");

    fwt53(im53, 1);
    fwt97(im97, 1);
    printf("level 1 Wavelets coefficients:\n");
    for (i = 0; i < n; i++) printf("wc[%d] = %f ; %f\n", i, im53[i], im97[i]);
    printf("\n");

    /*fwt97(im, 2);
    printf("level 2 Wavelets decomposition:\n");
    for (i = 0; i < n; i++) printf("wc[%d]=%f\n", i, im[i]);
    printf("\n");*/

    // Do the inverse 9/7 transform
    iwt53(im53, 1);
    iwt53(im53, 0);
    iwt97(im97, 1);
    iwt97(im97, 0);

    // Prints the reconstructed signal 
    printf("Reconstructed signal:\n");
    for (i = 0; i < n; i++) printf("im[%d] = %f , %f\n", i, im53[i], im97[i]);
}
