// exercise.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <vector>

int main()
{
	std::vector<double> image{ 7, 10, 8, 6, 4, 1, 3, 7 };
	size_t len = image.size();
	std::vector<std::vector<double>> decomp;
	std::vector<double> a(image);
	int level = 1;
	while (len > 1)
	{
		std::vector<double> d(len);
		len = len / 2;
		for (int i = 0; i < len; ++i)
		{
			double even = a[2 * i];
			double odd = a[2 * i + 1];
			d[i] = (even + odd) / 2;
			d[i + len] = (even - odd) / 2;
		}
		decomp.emplace_back(d);
		std::cout << "level " << level << "\n";
		for (int i = 0; i < 2 * len; ++i)
		{
			a[i] = d[i];
			std::cout << i << ": " << d[i] << "\n";
		}
		++level;
	}
	return 0;
}
