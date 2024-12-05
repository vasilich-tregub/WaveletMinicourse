// Le Gall 5/3 1D code.cpp : Defines the entry point for the application.
//

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
	// successive convolutions with {-1./4, 1., -1./4} for even pixels
	// and {1./2, 1., 1./2} for even pixels
	// for im[n] result is -im[n-2]/8 + im[n-1]/8 + 6*im[n]/8 + im[n+1]/8 - im[n+2]/8
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
	// for im[n] result is -im[n-2]/8 + im[n-1]/4 + 6*im[n]/8 + im[n+1]/4 - im[n+2]/8
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
	std::vector<int32_t> im{ 7, 10, 8, 6, 4, 1, 3, 7 };
	size_t len = im.size();

	for (int i = 0; i < len; ++i)
		im[i] *= 256 * 16;

	dwt_forward(im, 0);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";
	dwt_forward(im, 1);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";
	dwt_forward(im, 2);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";
	std::cout << "Inverse:\n";
	dwt_inverse(im, 2);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";
	dwt_inverse(im, 1);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";
	dwt_inverse(im, 0);
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] / 256.0 / 16 << "\n";


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
