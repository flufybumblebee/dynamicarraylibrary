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
	NeuralNetwork nn = { 2,2,1 };

	TrainningData td[4];
	td[0].inputs	= { 0.0f,0.0f};
	td[0].answers	= { 0.0f };
	td[1].inputs	= { 0.0f,1.0f };
	td[1].answers	= { 1.0f };
	td[2].inputs	= { 1.0f,0.0f };
	td[2].answers	= { 1.0f };
	td[3].inputs	= { 1.0f,1.0f };
	td[3].answers	= { 0.0f };

	for (size_t i = 0; i < 50000; i++)
	{	
		auto data = td[RandomInt(0, 3)];
		nn.Train(data.inputs, data.answers);		
	}

	for (int i = 0; i < 4; i++)
	{
		auto outputs = nn.Predict(td[i].inputs);
		std::cout << outputs.at(0) << std::endl;
	}

	return 0;
}
