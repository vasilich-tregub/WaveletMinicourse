// HAAR-interleaved.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <vector>

int main()
{
	std::vector<int> im{ 7, 10, 8, 6, 4, 1, 3, 7 };
	std::cout << "Source:" << "\n";
	for (int i = 0; i < im.size(); ++i)
	{
		std::cout << '#' << i << ":" << im[i] << "; ";
	}
	std::cout << "\n";
	size_t len = im.size();
	std::vector<std::vector<int>> decomp;
	std::cout << "Direct biorthogonal Haar transform:" << "\n";
	std::vector<int> a(im);
	int level = 1;
	while (len > 1)
	{
		std::vector<int> d(len);
		len = len / 2;
		for (int i = 0; i < len; ++i)
		{
			int even = a[2 * i];
			int odd = a[2 * i + 1];
			d[i] = (even + odd);			// Note that values (even + odd) and (even - odd) have identical parity
			d[i + len] = (even - odd);		// d[i] (approx) and d[i + len] (coeff) are either both even or both odd
		}
		decomp.emplace_back(d);
		std::cout << "level " << level << ":\n";
		for (int i = 0; i < 2 * len; ++i)
		{
			a[i] = d[i];
			std::cout << '#' << i << ":" << d[i] << "; ";
		}
		std::cout << "\n";
		++level;
	}
	std::cout << "\n";
	std::cout << "Decomposition:" << "\n";
	for (int i = 0; i < im.size(); ++i)
	{
		std::cout << '#' << i << ":" << a[i] << "; ";
	}
	std::cout << "\n" << "\n";
	std::cout << "Inverse biorthogonal Haar transform:" << "\n";
	std::vector<int> recovered(im.size());
	--level;
	while (level > 0)
	{
		len = decomp[--level].size();
		std::vector<int> recovered_coeffs(len);
		len = len / 2;
		for (int i = 0; i < len; ++i)
		{
			int approx = a[i];
			int coeff = a[i + len];
			recovered_coeffs[2 * i] = (approx + coeff) / 2;			// see lines 29, 30. 'approx' and 'coeff' values
			recovered_coeffs[2 * i + 1] = (approx - coeff) / 2;		// of the same decomp level have identical parity.
		}
		std::cout << "level " << level << "(recovered values):\n";
		for (int i = 0; i < 2 * len; ++i)
		{
			a[i] = recovered_coeffs[i];
			std::cout << '#' << i << ":" << recovered_coeffs[i] << "; ";
		}
		std::cout << "\n";
	}
	return 0;
}
