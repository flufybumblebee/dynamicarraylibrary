#pragma once

#include "Matrix.h"
#include "MathLibrary.h"

class NeuralNetwork
{
public:
	int input_nodes;
	int hidden_nodes;
	int output_nodes;
	Matrix weight_ih;
	Matrix weight_ho;
	Matrix bias_h;
	Matrix bias_o;

public:
	NeuralNetwork(
		const int& input_nodes,
		const int& hidden_nodes,
		const int& output_nodes)
		:
		input_nodes(input_nodes),
		hidden_nodes(hidden_nodes),
		output_nodes(output_nodes),
		weight_ih( Matrix( hidden_nodes, input_nodes )),
		weight_ho( Matrix( output_nodes, hidden_nodes)),
		bias_h( Matrix(hidden_nodes, 1)),
		bias_o(Matrix(output_nodes, 1))
	{
		weight_ih.Randomise(-1.0f, 1.0f);
		weight_ho.Randomise(-1.0f, 1.0f);
		bias_h.Randomise(-1.0f, 1.0f);
		bias_o.Randomise(-1.0f, 1.0f);
	}

public:
	std::vector<float> FeedForward(std::vector<float> input_array)
	{
		auto inputs = Matrix::FromArray(input_array);
		auto hidden = weight_ih * inputs;
		hidden += bias_h;
		hidden.Map(Sigmoid);

		auto outputs = Matrix::Multiply(weight_ho, hidden);
		outputs += bias_o;
		outputs.Map(Sigmoid);

		return outputs.ToArray();
	}

	void Train()
	{

	}
};
