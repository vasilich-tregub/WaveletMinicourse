// cdf 9/7.cpp : Defines the entry point for the application.
// based off the code snippet from https://gist.github.com/i-e-b/bb72fed460418f7c7ccb221d4b1da2b1
// credited, in turn, to '2006 - Gregoire Pau - gregoire.pau@ebi.ac.uk'
// I added a decomposition level to the parameter list and omitted decomposition packing/unpacking

/**
 *  Fast discrete biorthogonal CDF 9/7 wavelet forward and inverse transform (lifting implementation)
 *  dwt97.c 
 *  2006 - Gregoire Pau - gregoire.pau@ebi.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>

/**
 *  fwt97 - Forward biorthogonal 9/7 wavelet transform (lifting implementation)
 *
 *  im is an input signal, which will be replaced by its output transform.
 */

void fwt97(std::vector<double>& im, const int level) {
    const int inc = (int)1 << level;
    int end = (int)im.size();
    assert(inc < end && "stepping outside source image");

    double a;
    int i;

    // Predict 1
    a = -1.586134342;
    for (i = inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Update 1
    a = -0.05298011854;
    im[0] += 2 * a * im[inc];
    for (i = 2 * inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Predict 2
    a = 0.8829110762;
    for (i = inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Update 2
    a = 0.4435068522;
    im[0] += 2 * a * im[inc];
    for (i = 2 * inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Scale
    a = 1 / 1.149604398;
    for (i = 0; i < end; i++) {
        if (i % (2 * inc)) im[i] *= a;
        else im[i] /= a;
    }

}

/**
 *  iwt97 - Inverse biorthogonal 9/7 wavelet transform
 *
 *  This is the inverse of fwt97 so that iwt97(fwt97(im,0),0)=im
 */

void iwt97(std::vector<double>& im, const int level) {
    const int inc = (int)1 << level;
    int end = (int)im.size();
    assert(inc < end && "stepping outside source image");

    double a;
    int i;

    // Undo scale
    a = 1.149604398;
    for (i = 0; i < end; i++) {
        if (i % (2 * inc)) im[i] *= a;
        else im[i] /= a;
    }

    // Undo update 2
    a = -0.4435068522;
    im[0] += 2 * a * im[inc];
    for (i = 2 * inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end) 
        im[i] += 2 * a * im[i - inc];

    // Undo predict 2
    a = -0.8829110762;
    for (i = inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Undo update 1
    a = 0.05298011854;
    im[0] += 2 * a * im[inc];
    for (i = 2 * inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];

    // Undo predict 1
    a = 1.586134342;
    for (i = inc; i < end - inc; i += 2 * inc) {
        im[i] += a * (im[i - inc] + im[i + inc]);
    }
    if (i < end)
        im[i] += 2 * a * im[i - inc];
}

int main() {
    const int n = 64;
    std::vector<double> im(n);
    int i;

    // Makes a polinomial signal
    for (i = 0; i < n; i++) im[i] = 5 + i + 0.4 * i * i;// -0.02 * i * i * i;

    // Prints original signal im
    printf("Original signal:\n");
    for (i = 0; i < n; i++) printf("im[%d]=%f\n", i, im[i]);
    printf("\n");

    // Do the forward 9/7 transform
    fwt97(im, 0);

    // Prints the wavelet coefficients
    /*printf("level 0 Wavelets coefficients:\n");
    for (i = 0; i < n; i++) printf("wc[%d]=%f\n", i, im[i]);
    printf("\n");*/

    fwt97(im, 1);
    /*printf("level 1 Wavelets coefficients:\n");
    for (i = 0; i < n; i++) printf("wc[%d]=%f\n", i, im[i]);
    printf("\n");*/

    fwt97(im, 2);
    printf("level 2 Wavelets decomposition:\n");
    for (i = 0; i < n; i++) printf("wc[%d]=%f\n", i, im[i]);
    printf("\n");

    // Do the inverse 9/7 transform
    iwt97(im, 2);
    iwt97(im, 1);
    iwt97(im, 0);

    // Prints the reconstructed signal 
    printf("Reconstructed signal:\n");
    for (i = 0; i < n; i++) printf("im[%d]=%f\n", i, im[i]);
}
