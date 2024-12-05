#pragma once
#include <vector>
#include <assert.h>

void iwt53(std::vector<double>& im, const int level)
{
	const int inc = (int)1 << level;
	int end = (int)im.size();
	assert(inc < end && "stepping outside source image");

	// low pass filter, {-1./4, 1./4, -1./4}
	int i = 0;
	im[i] -= (im[inc] + 1) / 2;
	i += 2 * inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] -= (im[i - inc] + im[i + inc] + 2) / 4;
	}
	if (i < end)
	{
		im[i] -= (im[i - inc] + 1) / 2;
	}

	// high pass filter, {-1./8, 1./8, 6./8, 1./8 -1./8}
	// successive convolutions with {-1./4, 1./4, -1./4} for even pixels
	// and {1./2, 1., 1./2} for even pixels
	// for im[n] result is -im[n-2]/8 + im[n-1]/8 + 6*im[n]/8 + im[n+1]/8 - im[n+2]/8
	i = inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] += (im[i - inc] + im[i + inc]) / 2;
	}
	if (i < end)
	{
		im[i] += im[i - inc];
	}

}

void fwt53(std::vector<double>& im, const int level)
{
	const int inc = (int)1 << level;
	int end = (int)im.size();
	assert(inc < end && "stepping outside source image");

	int i = inc;
	// high pass filter, {-1./2, 1., -1./2}
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] -= (im[i - inc] + im[i + inc]) / 2;
	}
	if (i < end)
	{
		im[i] -= im[i - inc];
	}

	i = 0;
	// low pass filter, 
	// successive convolutions with {-1./2, 1., -1./2} for odd pixels
	// and {1./4, 1., 1./4} for even pixels
	// for im[n] result in -im[n-2]/8 + im[n-1]/4 + 6*im[n]/8 + im[n+1]/4 - im[n+2]/8
	// i.e., {-1./8, 2./8, 6./8, 2./8, -1./8}
	im[i] += (im[inc] + 1) / 2;
	i += 2 * inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] += (im[i - inc] + im[i + inc] + 2) / 4;
	}
	if (i < end)
	{
		im[i] += (im[i - inc] + 1) / 2;
	}
}

