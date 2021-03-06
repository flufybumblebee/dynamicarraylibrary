#include <iostream>
#include "NeuralNetwork.h"

class TrainningData
{
public:
	std::vector<float> inputs;
	std::vector<float> answers;
};

int main()
{
	NeuralNetwork nn = { 3,5,2 };

	const int NUM = 8;
	TrainningData td[NUM];

	td[0].inputs	= { 0.0f,0.0f,0.0f};
	td[0].answers	= { 0.0f,0.0f };
	td[1].inputs	= { 0.0f,0.0f,1.0f };
	td[1].answers	= { 0.0f,1.0f };
	td[2].inputs	= { 0.0f,1.0f,0.0f };
	td[2].answers	= { 0.0f,1.0f };
	td[3].inputs	= { 0.0f,1.0f,1.0f };
	td[3].answers	= { 1.0f,0.0f };

	td[4].inputs	= { 1.0f,0.0f,0.0f };
	td[4].answers	= { 0.0f,1.0f };
	td[5].inputs	= { 1.0f,0.0f,1.0f };
	td[5].answers	= { 1.0f,0.0f };
	td[6].inputs	= { 1.0f,1.0f,0.0f };
	td[6].answers	= { 1.0f,0.0f };
	td[7].inputs	= { 1.0f,1.0f,1.0f };
	td[7].answers	= { 1.0f,1.0f };

	for (size_t i = 0; i < 10000; i++)
	{	
		auto data = td[RandomInt(0, NUM-1)];
		nn.Train(data.inputs, data.answers);		
	}
	std::cout << "Cin  B    A    Cout Sum" << std::endl;
	for (int i = 0; i < NUM; i++)
	{
		auto inputs = td[i].inputs;
		auto outputs = nn.Predict(td[i].inputs);
		std::cout << int(inputs.at(0) + 0.5f) << "    ";
		std::cout << int(inputs.at(1) + 0.5f) << "    ";
		std::cout << int(inputs.at(2) + 0.5f) << "    ";
		std::cout << int(outputs.at(0) + 0.5f) << "    ";
		std::cout << int(outputs.at(1) + 0.5f) << std::endl;
	}
	// BINARY FULL ADDER
	return 0;
}
