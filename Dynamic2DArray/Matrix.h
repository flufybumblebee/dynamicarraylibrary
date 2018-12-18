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
		assert(a.nCols == b.nRows);
		Matrix result( a.nRows, b.nCols );
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
		return result;

		/* 
			notes:
			always reading data from a as row matrix
			always reading data from b as column matrix
			always returns a row matrix
		*/
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
		return *this = temp;

		/*
			notes:
			always reading data from a as row matrix
			always reading data from b as column matrix
			always returns a row matrix
		*/
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

	Matrix Transpose() const
	{
		Matrix result(nCols, nRows);
		for (int i = 0; i < nRows; i++)
		{
			for (int j = 0; j < nCols; j++)
			{
				result.data[j][i] = data[i][j];
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
		return *this = temp;
	}
	static Matrix Transpose( const Matrix& m )
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

