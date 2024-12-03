#pragma once
#include <vector>
#include <assert.h>

/**
 *  fwt97 - Forward biorthogonal 9/7 wavelet transform (lifting implementation)
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
