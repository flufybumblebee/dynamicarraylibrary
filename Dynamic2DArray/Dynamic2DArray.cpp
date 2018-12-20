#include <iostream>
#include "Matrix.h"

float Double(float x)
{
	return x * 2;
}

float Triple(float x)
{
	return x * 3;
}

int main()
{
	Matrix a(3, 3);
	
	a.Randomise(0,2);
	a.Print();
	std::cout << std::endl;

	Matrix b(3, 3);

	b.Randomise(0, 2);
	b.Print();
	std::cout << std::endl;
	
	Matrix c = a * b;
	c.Print();
	std::cout << std::endl;
	
	return 0;
}
