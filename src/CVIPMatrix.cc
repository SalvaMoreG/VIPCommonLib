
#include "CVIPMatrix.h"

#include "VIPconstants.h"
#include "CVIPUtils.h"

#include <cassert>
#include <cmath>


ostream& operator<<(ostream& os, const CVIPMatrix& mat)
{
	for (int irow=0; irow < mat.GetNumberOfRows(); irow++)
	{
		for (int jcol=0; jcol < mat.GetNumberOfColumns(); jcol++)
		{
			os << mat.GetElement( irow, jcol ) << "  ";
		}
		os << endl;
	}

    return os;
}

// CVIPMatrix methods ---------------------------------

CVIPMatrix::CVIPMatrix(int in_nrows, int in_ncolumns)
	: m_data(NULL)
	, m_nrows(in_nrows)
	, m_ncolumns(in_ncolumns)
{
	InitMatrixElements();
}

CVIPMatrix::CVIPMatrix( int in_nrows, int in_ncolumns, const double& init_value )
	: m_data(NULL)
	, m_nrows(in_nrows)
	, m_ncolumns(in_ncolumns)
{
    InitMatrixElements(init_value);
}


// copy constructor
CVIPMatrix::CVIPMatrix(const CVIPMatrix& in_obj)
    : m_data(NULL)
    , m_nrows(in_obj.GetNumberOfRows())
    , m_ncolumns(in_obj.GetNumberOfColumns())
{
	InitMatrixElements();

	*this = in_obj;
}

// assign operator
CVIPMatrix&
CVIPMatrix::operator= (const CVIPMatrix& in_obj)
{
	assert (m_data != NULL);
	assert (m_nrows == in_obj.GetNumberOfRows());
	assert (m_ncolumns == in_obj.GetNumberOfColumns());

	if (this != &in_obj)
	{
		for (int irow=0; irow < m_nrows; irow++)
		{
			assert (m_data[irow] != NULL);
			for (int jcol=0; jcol < m_ncolumns; jcol++)
			{
				m_data[irow][jcol] = in_obj.GetElement(irow, jcol);
			}
		}
	}
	return *this;
}

// destructor
CVIPMatrix::~CVIPMatrix()
{
	for (int irow=0; irow < m_nrows; irow++)
		delete [] m_data[irow];
	delete [] m_data;
}

bool
CVIPMatrix::operator==(const CVIPMatrix &in_obj) const
{
	assert (m_data != NULL);
	
	if (this == &in_obj) 
		return true;

	assert (m_nrows == in_obj.GetNumberOfRows());
	assert (m_ncolumns == in_obj.GetNumberOfColumns());

	for (int irow=0; irow < m_nrows; irow++)
	{
		assert (m_data[irow] != NULL);
		for (int jcol=0; jcol < m_ncolumns; jcol++)
		{
			if ( !doubleEquals( m_data[irow][jcol], in_obj.GetElement(irow, jcol), 0.0001) )
			{
				return false;
			}
		}
	}
    return true;
}

bool
CVIPMatrix::operator!=(const CVIPMatrix &in_obj) const
{
    return !(*this == in_obj);
}

const double&
CVIPMatrix::GetElement( int irow, int jcol ) const
{
	assert (m_data != NULL);
	assert (irow >= 0 && irow < m_nrows); 
	assert (jcol >= 0 && jcol < m_ncolumns); 
	return m_data[irow][jcol];
}

void
CVIPMatrix::SetElement( int irow, int jcol, const double& in_value )
{
	assert (m_data != NULL);
	assert (irow >= 0 && irow < m_nrows); 
	assert (jcol >= 0 && jcol < m_ncolumns); 
	m_data[irow][jcol] = in_value;
}

// operators
CVIPVector
CVIPMatrix::operator*(const CVIPVector& in_vec) const
{
	assert (m_data != NULL);
	assert (m_ncolumns == in_vec.GetSize());

	CVIPVector result(m_nrows);
	for (int irow = 0; irow < m_nrows; irow++)
	{
		result[irow] = 0.0;
		for (int jcol = 0; jcol < m_ncolumns; jcol++)
		{
			result[irow] += m_data[irow][jcol] * in_vec[jcol];
		}
	}

	return result;
}

CVIPMatrix
CVIPMatrix::operator+(const CVIPMatrix& in_m) const
{
	assert (m_data != NULL);
    assert (m_nrows == in_m.GetNumberOfRows());
	assert (m_ncolumns == in_m.GetNumberOfColumns());

	CVIPMatrix result(m_nrows, m_ncolumns);
	for (int irow = 0; irow < m_nrows; irow++)
	{
		for (int jcol = 0; jcol < m_ncolumns; jcol++)
        {
            double A = m_data[irow][jcol];
            double B = in_m.GetElement(irow, jcol);
            result.SetElement(irow, jcol, A+B);
		}
	}
	return result;
}

CVIPMatrix
CVIPMatrix::operator-(const CVIPMatrix& in_m) const
{
	assert (m_data != NULL);
    assert (m_nrows == in_m.GetNumberOfRows());
	assert (m_ncolumns == in_m.GetNumberOfColumns());

	CVIPMatrix result(m_nrows, m_ncolumns);
	for (int irow = 0; irow < m_nrows; irow++)
	{
		for (int jcol = 0; jcol < m_ncolumns; jcol++)
        {
            double A = m_data[irow][jcol];
            double B = in_m.GetElement(irow, jcol);
            result.SetElement(irow, jcol, A-B);
		}
	}
	return result;
}

void
CVIPMatrix::InitMatrixElements( const double& init_value )
{
	m_data = new double* [m_nrows];
	for (int irow=0; irow < m_nrows; irow++)
		m_data[irow] = new double [m_ncolumns];

	for (int irow=0; irow < m_nrows; irow++)
		for (int jcol=0; jcol < m_ncolumns; jcol++)
			m_data[irow][jcol] = init_value;
}




