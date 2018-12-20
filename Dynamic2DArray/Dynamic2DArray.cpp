#include <iostream>
#include "NeuralNetwork.h"

int main()
{
	NeuralNetwork nn = { 2,2,1 };

	std::vector<float> input = { 1.0f, 0.0f };

	auto output = nn.FeedForward(input);

	for (size_t i = 0; i < output.size(); i++)
	{
		std::cout << output[i] << std::endl;
	}
	
	return 0;
}
