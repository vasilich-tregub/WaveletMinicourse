// Le Gall 5/3 1D code.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <iomanip>
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
		im[i - inc] += (im[i - 2 * inc] + im[i]) >> 1;
	}
	if (i < end)
	{
		im[i] -= (im[i - inc] + 1) >> 1;
		im[i - inc] += (im[i - 2 * inc] + im[i]) >> 1;
	}
	else if (i - inc < end)
	{
		im[i - inc] += im[i - 2 * inc];
	}
	
	// high pass filter is included in the loop above (i - inc: lag-lead/lead-lag?),
	// successive convolutions with {-1./4, 1., -1./4} for even pixels
	// and {1./2, 1., 1./2} for odd pixels
	// for im[n] result is -im[n-2]/8 + im[n-1]/8 + 6*im[n]/8 + im[n+1]/8 - im[n+2]/8
	// i.e. {-1./8, 1./8, 6./8, 1./8 -1./8}
}

void dwt_forward(std::vector<int32_t>& im, const int level)
{
	const int inc = (int)1 << level;
	int end = (int)im.size();
	assert(inc < end && "stepping outside source image");

	int i = inc;
	// high pass filter, {-1./2, 1., -1./2}
	if (i >= end - inc)
	{
		im[i] -= im[i - inc];
		im[i - inc] += (im[i] + 1) >> 1;
		return;
	}
	im[i] -= (im[i - inc] + im[i + inc]) >> 1;
	im[i - inc] += (im[i] + 1) >> 1;
	i += 2 * inc;
	for (; i < end - 1 * inc; i += 2 * inc)
	{
		im[i] -= (im[i - inc] + im[i + inc]) >> 1; 
		im[i - inc] += (im[i - 2 * inc] + im[i] + 2) >> 2; 
	}
	if (i < end)
	{
		im[i] -= im[i - inc];
		im[i - inc] += (im[i - 2 * inc] + im[i] + 2) >> 2; 
	}
	else if (i - inc < end)
	{
		im[i - inc] += (im[i - 2 * inc] + 1) >> 1;
	}

	// low pass filter included in the loop above (i - inc: lead-lag/lag-lead?), 
	// successive convolutions with {-1./2, 1., -1./2} for odd pixels
	// and {1./4, 1., 1./4} for even pixels
	// for im[n] result is -im[n-2]/8 + im[n-1]/4 + 6*im[n]/8 + im[n+1]/4 - im[n+2]/8
	// i.e., {-1./8, 2./8, 6./8, 2./8, -1./8}
}


int main()
{
	std::vector<int32_t> im{ 7, 10, 8, 6, 4, 1, 3, 7, 15, 10, 8, 6, 4, 1, 3 };
	size_t len = im.size();

	for (int i = 0; i < len; ++i)
		std::cout << std::setw(3) << im[i] << "; ";
	std::cout << std::endl;

	std::cout << "Forward:\n";
	int level = 0;
	for (int k = 1; k < len; ++level, k *= 2)
	{
		dwt_forward(im, level);
		for (int i = 0; i < len; ++i)
			std::cout << std::setw(3) << im[i] << "; ";
		std::cout << std::endl;
	}
	level--; // back to max level of forward dwt operations
	std::cout << "Inverse:\n";
	for (int k = 1; k < len; level--, k *= 2)
	{
		dwt_inverse(im, level);
		for (int i = 0; i < len; ++i)
			std::cout << std::setw(3) << im[i] << "; ";
		std::cout << std::endl;
	}

	return 0;
}
