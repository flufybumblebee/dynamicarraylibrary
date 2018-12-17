#pragma once

#include <cassert>
#include "MathLibrary.h"

class Matrix
{
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
	Matrix& Randomise(const float& min, const float& max)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = RandomFloat(min,max);
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
				data[i][j] = float( RandomInt(min, max) );
			}
		}
		return *this;
	}
	
	// static functions (doesn't alter the matrix passed in)
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
	static Matrix Multiply(const Matrix& in, const float& num)
	{
		Matrix out = in;
		for (int i = 0; i < out.nRows; i++)
		{
			for (int j = 0; j < out.nCols; j++)
			{
				out.data[i][j] *= num;
			}
		}
		return out;
	}
	static Matrix DotProduct(const Matrix& a, const Matrix& b)
	{
		assert(a.nRows == b.nCols);
		Matrix temp = a;
		for (int j = 0; j < temp.nRows; j++)
		{
			for (int k = 0; k < temp.nCols; k++)
			{
				float sum = 0.0f;
				for (int i = 0; i < b.nRows; i++)
				{
					sum += a.data[j][i] * b.data[i][k];
				}
				temp.data[j][k] = sum;
			}
		}
		return temp;
	}
	// non-static functions (does alter the matrix)
	Matrix& Add( const float& num)
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
	Matrix& DotProduct(const Matrix& m)
	{
		assert(nRows == m.nCols);

		Matrix temp = *this;
		for (int j = 0; j < nRows; j++)
		{
			for (int k = 0; k < nCols; k++)
			{
				float sum = 0.0f;
				for (int i = 0; i < m.nCols; i++)
				{
					sum += data[j][i] * m.data[i][k];
				}
				temp.data[j][k] = sum;
			}
		}
		return *this = temp;
	}
	// operator overloading
	Matrix& operator = (const Matrix& in)
	{
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				data[i][j] = in.data[i][j];
			}
		}
		return *this;
	}
	Matrix operator + ( const float num ) const
	{
		Matrix out = *this;
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				out.data[i][j] += num;
			}
		}
		return out;
	}
	Matrix operator * (const float num) const
	{
		Matrix out = *this;
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				out.data[i][j] *= num;
			}
		}
		return out;
	}
public:
	// testing functions:
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
public:
	const int nRows;
	const int nCols;

	float** data = nullptr;
};

