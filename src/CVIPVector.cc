
#include "CVIPVector.h"

#include "CVIPUtils.h"

#include <cmath>
#include <cassert>


ostream& operator<<(ostream& os, const CVIPVector& vec)
{
	for (int index = 0; index < vec.GetSize(); index++)
	{
    	os << vec[index] << ' ' ;
	}
    return os;
}

// ===================================================================

CVIPVector::CVIPVector(int nvoxels)
	: m_nvoxels(nvoxels)
{
	m_data = new double [m_nvoxels];
	for (int i = 0; i < m_nvoxels; i++)
		m_data[i] = 0.0;
}

CVIPVector::CVIPVector( int nvoxels, const double& in_initvalue )
	: m_nvoxels(nvoxels)
{
	m_data = new double [m_nvoxels];
	for (int i = 0; i < m_nvoxels; i++)
		m_data[i] = in_initvalue;
}

CVIPVector::~CVIPVector()
{
	delete [] m_data;
}

// copy constructor
CVIPVector::CVIPVector(const CVIPVector& in_obj)
	: m_nvoxels(in_obj.GetSize())
{
	m_data = new double [m_nvoxels];
	for (int i = 0; i < m_nvoxels; i++)
		m_data[i] = 0.0;

	*this = in_obj;
}

// assignment operator
CVIPVector&
CVIPVector::operator= (const CVIPVector& in_obj)
{
	assert (m_nvoxels == in_obj.GetSize());
	if (this != &in_obj)
	{
		for (int index = 0; index < m_nvoxels; index++)
		{
			m_data[index] = in_obj[index];
		}
	}
	return *this;
}

double&
CVIPVector::operator[](int index)
{ 
	assert(index >= 0 && index < m_nvoxels); 
	return m_data[index]; 
}

const double&
CVIPVector::operator[](int index) const
{ 
	assert(index >= 0 && index < m_nvoxels); 
	return m_data[index]; 
}

bool
CVIPVector::operator==(const CVIPVector &in_obj) const
{
	if (m_nvoxels != in_obj.GetSize()) 
		return false;
	for (int index = 0; index < m_nvoxels; index++)
	{
		if (!doubleEquals( m_data[index], in_obj[index], 0.0001 ))
		{
			return false;
		}
	}
	return true;
}

bool
CVIPVector::operator!=(const CVIPVector &in_obj) const
{
    return !(*this == in_obj);
}

// operators
// + (plus) operator
CVIPVector
CVIPVector::operator+ (const CVIPVector& in_vec) const
{
	assert (m_nvoxels == in_vec.GetSize());
	CVIPVector result(m_nvoxels);
	for (int index = 0; index < m_nvoxels; index++)
	{
		result[index] = m_data[index] + in_vec[index];
	}
	return result;
}

// - (minus) operator
CVIPVector
CVIPVector::operator- (const CVIPVector& in_vec) const
{
	assert (m_nvoxels == in_vec.GetSize());
	CVIPVector result(m_nvoxels);
	for (int index = 0; index < m_nvoxels; index++)
	{
		result[index] = m_data[index] - in_vec[index];
	}
	return result;
}

// scalar product
double
CVIPVector::operator*(const CVIPVector& in_vec) const
{
	assert (m_nvoxels == in_vec.GetSize());
	double result;
	for (int index = 0; index < m_nvoxels; index++)
	{
		result += m_data[index] + in_vec[index];
	}
	return result;
}

CVIPVector
CVIPVector::operator*(double const& in_coeff) const
{
	CVIPVector result(m_nvoxels);
	for (int index = 0; index < m_nvoxels; index++)
	{
		result[index] = m_data[index] * in_coeff;
	}
	return result;
}












