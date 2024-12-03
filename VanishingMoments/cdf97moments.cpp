// cdf 9/7.cpp : Defines the entry point for the application.
// this exercise examines vanishing moments of CDF 9/7 wavelet
// the wavelet code is based off 
// the code snippet from https://gist.github.com/i-e-b/bb72fed460418f7c7ccb221d4b1da2b1
// credited, in turn, to '2006 - Gregoire Pau - gregoire.pau@ebi.ac.uk'
// I added a decomposition level to the parameter list, omitted decomposition packing/unpacking
// and corrected the scaling step which should be a 'leapfrog' summation

/**
 *  Fast discrete biorthogonal CDF 9/7 wavelet forward and inverse transform (lifting implementation)
 *  dwt97.c 
 *  2006 - Gregoire Pau - gregoire.pau@ebi.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <iostream>

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
    for (i = 0; i < end; i += inc) {
        if (i % (2 * inc)) im[i] /= a;
        else im[i] *= a;
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
    for (i = 0; i < end; i += inc) {
        if (i % (2 * inc)) im[i] /= a;
        else im[i] *= a;
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
    int32_t len = 64;
    int32_t maxval = 1024;
    std::vector<double> cst(len, maxval);
    std::vector<double> lin(len);
    for (int i = 0; i < len; ++i)
        lin[i] = i * maxval / len;
    std::vector<double> sqt(len);
    for (int i = 0; i < len; ++i)
        sqt[i] = (i - len / 2) * (i - len / 2) * maxval * 4 / len / len;
    std::vector<double> cub(len);
    for (int i = 0; i < len; ++i)
        cub[i] = (i - len / 2) * (i - len / 2) * (i - len / 2) * maxval * 8 / len / len / len;
    for (int i = 0; i < len; ++i)
    {
        cst[i] *= 256 * 16;
        lin[i] *= 256 * 16;
        sqt[i] *= 256 * 16;
        cub[i] *= 256 * 16;
    }
    fwt97(cst, 0);
    fwt97(lin, 0);
    fwt97(sqt, 0);
    fwt97(cub, 0);
    fwt97(sqt, 1);
    fwt97(cub, 1);
    for (int i = 0; i < len; ++i)
    {
        cst[i] /= 256 * 16;
        lin[i] /= 256 * 16;
        sqt[i] /= 256 * 16;
        cub[i] /= 256 * 16;
    }
    std::cout << "constant function details:\n";
    for (int i = 1; i < len; ++++i)
    {
        std::cout << cst[i] << "\n";
    }
    std::cout << "\nlinear function details:\n";
    for (int i = 1; i < len; ++++i)
    {
        std::cout << lin[i] << "\n";
    }
    std::cout << "\nsquare fwt:\n";
    for (int i = 0; i < len; ++i)
    {
        std::cout << sqt[i] << "\n";
    }
    std::cout << "\ncubic fwt:\n";
    for (int i = 0; i < len; ++i)
    {
        std::cout << cub[i] << "\n";
    }
    std::cout << std::endl;

    return 0;
}
