#pragma once

#include <random>

inline double RandomDouble(double min, double max)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> dist(min, max);

	return dist(rng);
}
inline float RandomFloat(float min, float max)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<float> dist(min, max);

	return dist(rng);
}
inline int RandomInt(int min, int max)
{
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> dist(min, max);

	return dist(rng);
}
inline int Sign(float num_in)
{
	if (num_in >= 0.0f)	return  1;
	else				return -1;
}

inline float Sigmoid(float x)
{
	return 1.0f / (1.0f + std::expf(-x));
}

inline float DSigmoid(float x)
{
	// actual sigmoid derivitive
	// return Sigmoid(x) * (1 - Sigmoid(x));
	// fake function because sigmoid is already implemented
	return x * (1.0f - x);
}