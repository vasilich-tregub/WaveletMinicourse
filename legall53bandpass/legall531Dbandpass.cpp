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
	// and {1./2, 1., 1./2} for even pixels
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
	std::vector<int32_t> im{ 7, 10, 8, 6, 4, 1, 3, 7, 15 };
	size_t len = im.size();

	for (int i = 0; i < len; ++i)
		im[i] *= 256 * 16;

	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	dwt_forward(im, 0);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	dwt_forward(im, 1);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	dwt_forward(im, 2);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	std::cout << "Inverse:\n";
	dwt_inverse(im, 2);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	dwt_inverse(im, 1);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;
	dwt_inverse(im, 0);
	for (int i = 0; i < len; ++i)
		std::cout << std::setw(6) << im[i] / 256.0 / 16 << "; ";
	std::cout << std::endl;


	/*std::vector<double> decomp(len);
	for (int d = 1; d <= len / 2; d *= 2)
	{
		for (int i = 0; i < len / 2 / d; ++i)
		{
			decomp[i] = im[d * i];
			decomp[i + len / 2 / d] = im[2 * d * i + d];
		}
	}
	std::cout << "Deinterleaved vector:\n";
	for (int i = 0; i < len; ++i)
	{
		std::cout << i << ": " << decomp[i] / 256 / 16 << "\n";
	}*/
	return 0;
}
