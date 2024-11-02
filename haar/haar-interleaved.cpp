// HAAR-interleaved.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <vector>

int main()
{
	std::vector<double> im{ 7, 10, 8, 6, 4, 1, 3, 7 };
	size_t len = im.size();
	int k = 1;
	while (k < len)
	{
		for (int i = 0; i < len; i += 2 * k)
		{
			double even = im[i];
			double odd = im[i + k];
			im[i] = (even + odd) / 2;
			im[i + k] = (even - odd) / 2;
		}
		k *= 2;
	}
	for (int i = 0; i < len; ++i)
		std::cout << i << ": " << im[i] << "\n";

/*
* "interleaved" decomposition created with above iterations:
* [5.75], -1.5,  0.75,  1,     2,  1.5, -1.25, -2
*
* "manual" decomposition to [approx], Hn, ... H1
* n = log2(len)
* [7,   10,      8,     6,     4,   1,    3,   7]  // image
* [8.5,  7,     2.5,    5],  -1.5,  1,   1.5, -2   // A       H1
* [7.75, 3.75], 0.75, -1.25, -1.5,  1,   1.5, -2   // A    H2 H1
* [5.75],   2,  0.75, -1.25, -1.5,  1,   1.5, -2   // A H2 H1 H3
* 
* compare to programmatic "de-interleave"
*/
	std::vector<double> decomp(len);
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
		std::cout << i << ": " << decomp[i] << "\n";
	}
	return 0;
}
