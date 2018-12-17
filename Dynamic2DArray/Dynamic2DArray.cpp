#include <iostream>
#include "Matrix.h"

int main()
{
	Matrix a(3, 2);
	Matrix b(2, 3);
	
	a.Randomise(1,3);
	a.Print();
	std::cout << std::endl;
	b.Randomise(1,3);
	b.Print();
	std::cout << std::endl;
	Matrix c = Matrix::Add(a, 3);
	c.Print();
	std::cout << std::endl;
	c.Add(3);
	std::cout << std::endl;
	c.Print();
	std::cout << std::endl;
	Matrix d = a + 5;
	std::cout << std::endl;
	d.Print();

	return 0;
}
