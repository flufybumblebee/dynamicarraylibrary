#pragma once

#include <cassert>
#include <memory>
#include <array>
#include "MathLibrary.h"

class Matrix
{
public:
	int nRows;
	int nCols;

	float** data = nullptr;
public:
	// constructors and destructors
	Matrix() = default;
	Matrix(const Matrix& m)
		:
		nRows(m.nRows),
		nCols(m.nCols)
	{
		data = new float*[nRows];

		for (int i = 0; i < nRows; i++)
		{
			data[i] = new float[nCols];
		}

		// copy data
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = m.data[i][j];
			}
		}
	}
	Matrix(const int nRows, const int nCols)
		:
		nRows(nRows),
		nCols(nCols)
	{
		data = new float*[nRows];

		for (int i = 0; i < nRows; i++)
		{
			data[i] = new float[nCols];
		}

		// set all to 0
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = 0.0f;
			}
		}
	}
	~Matrix()
	{
		for (int i = 0; i < nRows; i++)
		{
			delete[] data[i];
		}
		delete[] data;
		data = nullptr;
	}

	// functions
	static Matrix Transpose(const Matrix& m)
	{
		Matrix result(m.nCols, m.nRows);
		for (int i = 0; i < m.nRows; i++)
		{
			for (int j = 0; j < m.nCols; j++)
			{
				result.data[j][i] = m.data[i][j];
			}
		}
		return result;
	}
	Matrix& Transpose()
	{
		Matrix temp(nCols, nRows);
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				temp.data[j][i] = data[i][j];
			}
		}
		return *this = std::move(temp);
	}
	Matrix& Randomise(const float& min, const float& max)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = RandomFloat(min, max);
			}
		}
		return *this;
	}
	Matrix& Randomise(const int& min, const int& max)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = float(RandomInt(min, max));
			}
		}
		return *this;
	}
	static Matrix Map(float(*pFunc)(float), const Matrix& in)
	{
		Matrix out = in;
		for (int i = 0; i < out.nRows; i++)
		{
			for (int j = 0; j < out.nCols; j++)
			{
				auto var = out.data[i][j];
				out.data[i][j] = pFunc(var);
			}
		}
		return out;
	}
	Matrix& Map(float(*pFunc)(float))
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				auto var = data[i][j];
				data[i][j] = pFunc(var);
			}
		}
		return *this;
	}
	static Matrix FromArray(std::vector<float> input_array)
	{
		Matrix output = Matrix(input_array.size(), 1);

		for (size_t i = 0; i < input_array.size(); i++)
		{
			output.data[i][0] = input_array[i];
		}
		return std::move(output);
	}
	std::vector<float> ToArray() const
	{
		std::vector<float> array_out;
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				array_out.push_back(data[i][j]);
			}
		}
		return std::move(array_out);
	}

	// static functions (the matrix or matrices passed in aren't altered )
	static Matrix Add(const Matrix& in, const float& num)
	{
		Matrix out = in;
		for (int i = 0; i < out.nRows; i++)
		{
			for (int j = 0; j < out.nCols; j++)
			{
				out.data[i][j] += num;
			}
		}
		return out;
	}
	static Matrix Add(const Matrix& a, const Matrix& b)
	{
		assert(a.nRows == b.nRows && a.nCols == b.nCols);

		Matrix result(a.nRows, b.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = a.data[i][j] + b.data[i][j];
			}
		}
		return result;
	}
	static Matrix Subtract(const Matrix& m, const float& num)
	{
		Matrix result(m.nRows, m.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = m.data[i][j] - num;
			}
		}
		return std::move(result);
	}
	static Matrix Subtract(const Matrix& a, const Matrix& b)
	{
		assert(a.nRows == b.nRows && a.nCols == b.nCols);

		Matrix result(a.nRows,b.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = a.data[i][j] - b.data[i][j];
			}
		}
		return std::move(result);
	}
	static Matrix Multiply(const Matrix& m, const float& num)
	{
		Matrix result(m.nRows,m.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = m.data[i][j] * num;
			}
		}
		return std::move(result);
	}
	static Matrix Multiply(const Matrix& a, const Matrix& b)
	{
		assert(a.nCols == b.nRows);
		Matrix result(a.nRows, b.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				float sum = 0.0f;
				for (int k = 0; k < a.nCols; k++)
				{
					sum += a.data[i][k] * b.data[k][j];
				}
				result.data[i][j] = sum;
			}
		}
		return std::move(result);

		/*
			notes:
			always reading data from a as row matrix
			always reading data from b as column matrix
			always returns a row matrix
		*/
	}
	static Matrix Hadamard(const Matrix& a, const Matrix& b)
	{
		assert(a.nRows == b.nRows && a.nCols == b.nCols);
		Matrix result( a.nRows, b.nCols );
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = a.data[i][j] * b.data[i][j];
			}
		}
		return std::move(result);
	}
	
	// non-static functions (alters this matrix)
	Matrix& Add(const float& num)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] += num;
			}
		}
		return *this;
	}
	Matrix& Add(const Matrix& m)
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] += m.data[i][j];
			}
		}
		return *this;
	}
	Matrix& Subtract(const float& num)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] -= num;
			}
		}
		return *this;
	}
	Matrix& Subtract(const Matrix& m)
	{
		assert(nRows == m.nRows && nCols == m.nCols);
		
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] -= m.data[i][j];
			}
		}
		return *this;
	}
	Matrix& Multiply(const float& num)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] *= num;
			}
		}
		return *this;
	}
	Matrix& Multiply(const Matrix& m)
	{
		assert(nCols == m.nRows);
		Matrix temp(nRows, m.nCols);
		for (int i = 0; i < temp.nRows; i++)
		{
			for (int j = 0; j < temp.nCols; j++)
			{
				float sum = 0.0f;
				for (int k = 0; k < nCols; k++)
				{
					sum += data[i][k] * m.data[k][j];
				}
				temp.data[i][j] = sum;
			}
		}
		return *this = std::move(temp);

		/*
			notes:
			always reading data from a as row matrix
			always reading data from b as column matrix
			always returns a row matrix
		*/
	}
	Matrix& Hadamard(const Matrix& m)
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] *= m.data[i][j];
			}
		}
		return *this;
	}
	
	// operator overloading
	Matrix&	operator = (const Matrix& m)
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = m.data[i][j];
			}
		}
		return *this;
	}
	Matrix	operator + (const float num) const
	{
		Matrix result(nRows, nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = data[i][j] + num;
			}
		}
		return std::move(result);
	}
	Matrix	operator + (const Matrix& m) const
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		Matrix result(nRows,m.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = data[i][j] + m.data[i][j];
			}
		}
		return std::move(result);
	}
	Matrix&	operator += (const Matrix& m)
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] += m.data[i][j];
			}
		}
		return *this;
	}
	Matrix	operator - (const float num) const
	{
		Matrix result(nRows, nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = data[i][j] - num;
			}
		}
		return std::move(result);
	}
	Matrix	operator - (const Matrix& m) const
	{
		assert(nRows == m.nRows && nCols == m.nCols);

		Matrix result(nRows,m.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = data[i][j] - m.data[i][j];
			}
		}
		return std::move(result);
	}
	Matrix&	operator -= (const Matrix& m)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] -= m.data[i][j];
			}
		}
		return *this;
	}
	Matrix	operator * (const float num) const
	{
		Matrix result(nRows, nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				result.data[i][j] = data[i][j] * num;
			}
		}
		return std::move(result);
	}
	Matrix	operator * (const Matrix& m) const
	{
		assert(nCols == m.nRows);
		Matrix result(nRows, m.nCols);
		for (int i = 0; i < result.nRows; i++)
		{
			for (int j = 0; j < result.nCols; j++)
			{
				float sum = 0.0f;
				for (int k = 0; k < nCols; k++)
				{
					sum += data[i][k] * m.data[k][j];
				}
				result.data[i][j] = sum;
			}
		}
		return std::move(result);

		/*
			notes:
			always reading data from a as row matrix
			always reading data from b as column matrix
			always returns a row matrix
		*/
	}

public:
	// console testing functions:
	static void Print(const Matrix& m)
	{
		for (int i = 0; i < m.nRows; i++)
		{
			for (int j = 0; j < m.nCols; j++)
			{
				std::cout << m.data[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	void Print() const
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				std::cout << data[i][j] << "\t";
			}
			std::cout << std::endl;
		}
	}
};

