// legall53moments.cpp : Defines the entry point for the application.
// this exercise examines vanishing moments of Le Gall 5/3 wavelet
// the wavelet code is based off 
// the dwt code (dwt.h, dwt.c) of ISO 21122 - 5 ed.2 reference software
// credited, in turn, to Fraunhofer IIS and intopix.com
// I added comments emphasizing the use of lifting scheme implementation
// and perfect (bit-precise) reconstruction with lifting scheme

#include <iostream>
#include <vector>
#include <assert.h>

void dwt_inverse(std::vector<int32_t>& im, const int level)
{
	const int inc = (int)1 << level;
	int end = (int)im.size();
	assert(inc < end && "stepping outside source image");

	// low pass filter, {-1./4, 1./4, -1./4}
	int i = 0;
	im[i] -= (im[inc] + 1) >> 1;
	i += 2 * inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] -= (im[i - inc] + im[i + inc] + 2) >> 2;
	}
	if (i < end)
	{
		im[i] -= (im[i - inc] + 1) >> 1;
	}
	
	// high pass filter, {-1./8, 1./8, 6./8, 1./8 -1./8}
	// successive convolutions with {-1./4, 1./4, -1./4} for even pixels
	// and {1./2, 1., 1./2} for even pixels
	// for im[n] result in -im[n-2]/8 + im[n-1]/8 + 6*im[n]/8 + im[n+1]/8 - im[n+2]/8
	i = inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] += (im[i - inc] + im[i + inc]) >> 1;
	}
	if (i < end)
	{
		im[i] += im[i - inc];
	}

}

void dwt_forward(std::vector<int32_t>& im, const int level)
{
	const int inc = (int)1 << level;
	int end = (int)im.size();
	assert(inc < end && "stepping outside source image");

	int i = inc;
	// high pass filter, {-1./2, 1., -1./2}
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] -= (im[i - inc] + im[i + inc]) >> 1; 
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
	im[i] += (im[inc] + 1) >> 1;
	i += 2 * inc;
	for (; i < end - inc; i += 2 * inc)
	{
		im[i] += (im[i - inc] + im[i + inc] + 2) >> 2; 
	}
	if (i < end)
	{
		im[i] += (im[i - inc] + 1) >> 1;
	}
}


int main()
{
	uint32_t len = 32;
	uint32_t maxval = 1024;
	std::vector<int32_t> cst(len, maxval);
	std::vector<int32_t> lin(len);
	for (int i = 0; i < len; ++i)
		lin[i] = i * maxval / len;
	std::vector<int32_t> sqt(len);
	for (int i = 0; i < len; ++i)
		sqt[i] = (i - len / 2) * (i - len / 2) * maxval * 4 / len / len;
	std::vector<int32_t> cub(len);
	for (int i = 0; i < len; ++i)
		cub[i] = (i - len / 2) * (i - len / 2) * (i - len / 2) * maxval * 8 / len / len / len;
	for (int i = 0; i < len; ++i)
	{
		cst[i] *= 256 * 16;
		lin[i] *= 256 * 16;
		sqt[i] *= 256 * 16;
		cub[i] *= 256 * 16;
	}
	dwt_forward(cst, 0);
	dwt_forward(lin, 0);
	dwt_forward(sqt, 0);
	dwt_forward(cub, 0);
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
		std::cout << cst[i] << " ";
	}
	std::cout << "\nlinear function details:\n";
	for (int i = 1; i < len; ++++i)
	{
		std::cout << lin[i] << " ";
	}
	std::cout << "\nsquare polynomial details:\n";
	for (int i = 1; i < len; ++++i)
	{
		std::cout << sqt[i] << " ";
	}
	std::cout << "\ncubic polynomial details:\n";
	for (int i = 1; i < len; ++++i)
	{
		std::cout << cub[i] << " ";
	}
	std::cout << std::endl;
	return 0;
}
