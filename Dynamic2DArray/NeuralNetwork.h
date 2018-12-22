#pragma once

#include "Matrix.h"
#include "MathLibrary.h"

class NeuralNetwork
{
public:
	int input_nodes;
	int hidden_nodes;
	int output_nodes;
	Matrix weights_ih;
	Matrix weights_ho;
	Matrix bias_h;
	Matrix bias_o;
	const float learning_rate = 0.2f;

public:
	NeuralNetwork(
		const int& input_nodes,
		const int& hidden_nodes,
		const int& output_nodes)
		:
		input_nodes(input_nodes),
		hidden_nodes(hidden_nodes),
		output_nodes(output_nodes),
		weights_ih( Matrix( hidden_nodes, input_nodes )),
		weights_ho( Matrix( output_nodes, hidden_nodes)),
		bias_h( Matrix(hidden_nodes, 1)),
		bias_o(Matrix(output_nodes, 1))
	{
		weights_ih.Randomise(-1.0f, 1.0f);
		weights_ho.Randomise(-1.0f, 1.0f);
		bias_h.Randomise(-1.0f, 1.0f);
		bias_o.Randomise(-1.0f, 1.0f);
	}

public:
	std::vector<float> Predict(std::vector<float> input_array)
	{
		// hidden weights
		auto inputs = Matrix::FromArray(input_array);
		auto hidden = weights_ih * inputs;
		hidden += bias_h;
		hidden.Map(Sigmoid);

		// output weights
		auto outputs = Matrix::Multiply(weights_ho, hidden);
		outputs += bias_o;
		outputs.Map(Sigmoid);

		return outputs.ToArray();
	}

	void Train(std::vector<float> inputs_array,std::vector<float> answers_array)
	{
		// hidden weights
		auto inputs = Matrix::FromArray(inputs_array);
		auto hidden = weights_ih * inputs;
		hidden += bias_h;
		hidden.Map(Sigmoid);

		// output weights
		auto outputs = weights_ho * hidden;
		outputs += bias_o;
		outputs.Map(Sigmoid);

		// calculate output layer errors
		auto answers = Matrix::FromArray(answers_array);
		auto output_errors = answers - outputs;
		
		// calculate gradient
		auto gradients = Matrix::Map(DSigmoid, outputs);
		gradients.Hadamard(output_errors);
		gradients.Multiply(learning_rate);

		// calculate deltas
		auto hidden_t = Matrix::Transpose(hidden);
		auto weight_ho_deltas = gradients * hidden_t;
		
		// adjust weights by deltas
		weights_ho.Add(weight_ho_deltas);

		// adjust the bias by its deltas (which is the gradients)
		bias_o.Add(gradients);

		// calculate hidden layer errors
		auto weight_ho_t = Matrix::Transpose(weights_ho);
		auto hidden_errors = weight_ho_t * output_errors;

		// calculate hidden gradient
		auto hidden_gradient = Matrix::Map(DSigmoid, hidden);
		hidden_gradient.Hadamard(hidden_errors);
		hidden_gradient.Multiply(learning_rate);

		// calculate input to hidden deltas
		auto inputs_t = Matrix::Transpose(inputs);
		auto weights_ih_deltas = hidden_gradient * inputs_t;

		// adjust weights by deltas
		weights_ih.Add(weights_ih_deltas);

		// adjust the bias by its deltas(which is the gradients)
		bias_h.Add(hidden_gradient);
	}
};
