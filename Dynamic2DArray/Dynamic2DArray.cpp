#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix a(1, 3);
	Matrix b(3, 3);
	
	a.Randomise(1,3);
	a.Print();
	std::cout << std::endl;
	b.Randomise(1,3);
	b.Print();
	std::cout << std::endl;
	Matrix c = Matrix::Transpose(b);
	c.Print();
	std::cout << std::endl;
	b.Print();

	return 0;
}
