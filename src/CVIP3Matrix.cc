
#include "CVIP3Matrix.h"
#include "VIPconstants.h"

#include <cassert>
#include <cmath>

using namespace std;

ostream& operator<<(ostream& os, const C3Matrix& mat)
{
	C3Vector tmp;
	mat.GetRow(0, tmp);
    os << tmp << endl;
	mat.GetRow(1, tmp);
    os << tmp << endl;
	mat.GetRow(2, tmp);
    os << tmp;

    return os;
}

// C3Matrix methods ---------------------------------

C3Matrix::C3Matrix()
	: m_a(NULL)
{
	InitMatrixElements();
}


C3Matrix::C3Matrix(const C3Vector& in_col1, const C3Vector& in_col2, const C3Vector& in_col3)
	: m_a(NULL)
{
	InitMatrixElements();

	Set(in_col1, in_col2, in_col3);
}

C3Matrix::C3Matrix( AXIS in_axis, const double& in_angle)
{
	InitMatrixElements();

	Set( in_axis, in_angle );
}

// copy constructor
C3Matrix::C3Matrix(const C3Matrix& in_obj)
{
	InitMatrixElements();

	*this = in_obj;
}


// assign operator
C3Matrix&
C3Matrix::operator= (const C3Matrix& in_obj)
{
	assert (m_a != NULL);
	if (this != &in_obj)
	{
		for (int icol=0; icol < 3; icol++)
		{
			for (int jrow=0; jrow < 3; jrow++)
			{
				m_a[icol][jrow] = in_obj.GetElement(icol, jrow);
			}
		}
	}
	return *this;
}

// destructor
C3Matrix::~C3Matrix()
{
	for (int i=0; i < 3; i++)
		delete [] m_a[i];
	delete [] m_a;
}


bool
C3Matrix::operator==(const C3Matrix &in_obj) const
{
	if (this == &in_obj) return 1;

	for (int icol=0; icol < 3; icol++)
	{
		for (int jrow=0; jrow < 3; jrow++)
		{
			if ( !doubleEquals( m_a[icol][jrow], in_obj.GetElement(icol, jrow), 0.0001) )
			{
				return 0;
			}
		}
	}
    return 1;
}

bool
C3Matrix::operator!=(const C3Matrix &in_obj) const
{
    return !(*this == in_obj);
}

void
C3Matrix::Set (const C3Vector& in_col1, const C3Vector& in_col2, const C3Vector& in_col3)
{
	assert (m_a != NULL);

	const C3Vector* vec;
	for (int icol=0; icol < 3; icol++)
	{
		if (icol==0) vec = &in_col1;
		else if (icol==1) vec = &in_col2;
		else if (icol==2) vec = &in_col3;

		m_a[icol][0] = vec->GetX();
		m_a[icol][1] = vec->GetY();
		m_a[icol][2] = vec->GetZ();
	}
}

void
C3Matrix::Set( AXIS in_axis, const double& in_angle)
{
	double sinAngle = sin(in_angle);
	double cosAngle = cos(in_angle);

	if (in_axis == Axis_X)
	{
		// cout << "X angle sin/cos angle: " << in_angle << " " << sinAngle << " " << cosAngle << endl;
		m_a[0][0] = 1.0;  m_a[1][0] = 0.0; 			m_a[2][0] = 0.0;
		m_a[0][1] = 0.0;  m_a[1][1] = cosAngle; 	m_a[2][1] = -sinAngle;
		m_a[0][2] = 0.0;  m_a[1][2] = sinAngle; 	m_a[2][2] = cosAngle;
	}
	else if (in_axis == Axis_Y)
	{
		// cout << "Y angle sin/cos angle: " << in_angle << " " << sinAngle << " " << cosAngle << endl;
		m_a[0][0] = cosAngle;  		m_a[1][0] = 0.0; 	m_a[2][0] = sinAngle;
		m_a[0][1] = 0.0;  			m_a[1][1] = 1; 		m_a[2][1] = 0.0;
		m_a[0][2] = -sinAngle;  	m_a[1][2] = 0.0; 	m_a[2][2] = cosAngle;
	}
	else if (in_axis == Axis_Z)
	{
		// cout << "Z angle sin/cos angle: " << in_angle << " " << sinAngle << " " << cosAngle << endl;
		m_a[0][0] = cosAngle;  	m_a[1][0] = -sinAngle;  m_a[2][0] = 0.0;
		m_a[0][1] = sinAngle;  	m_a[1][1] = cosAngle; 	m_a[2][1] = 0.0;
		m_a[0][2] = 0.0;  		m_a[1][2] = 0.0; 		m_a[2][2] = 1.0;
	}
	// cout << "MATRIX: " << *this << endl;
}

double
C3Matrix::GetElement(int col, int row) const
{
	assert (m_a != NULL);
	assert (*m_a != NULL);

	return m_a[col][row];
}

void
C3Matrix::SetElement(int col, int row, const double& value)
{
	assert (m_a != NULL);

	m_a[col][row] = value;
}

void
C3Matrix::SetUnity( const double& in_factor )
{
	for (int icol=0; icol < 3; icol++)
	{
		for (int jrow=0; jrow < 3; jrow++)
		{
			if (icol == jrow)
				m_a[icol][jrow] = in_factor;
			else
				m_a[icol][jrow] = 0.0;
		}
	}
}

void
C3Matrix::GetColumn(int in_col, C3Vector& io_vec) const
{
	assert (m_a != NULL);
	assert (in_col >= 0 && in_col < 3);

	io_vec.Set( m_a[in_col][0], m_a[in_col][1], m_a[in_col][2] );
}

void
C3Matrix::GetRow(int in_row, C3Vector& io_vec) const
{
	assert (m_a != NULL);
	assert (in_row >= 0 && in_row < 3);

	io_vec.Set( m_a[0][in_row], m_a[1][in_row], m_a[2][in_row] );
}

/*
double
C3Matrix::GetInverse() const
{
	// pffff... Not implemented...
	assert (false);
}
*/

bool		
C3Matrix::IsSymmetric() const
{
	// TODO: is there a smarter algorithm, faster than O(N^2)?
	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			if (j != i && !doubleEquals( m_a[i][j], m_a[j][i], 0.000001 ))
			{
				return false;
			}				
		}
	}
	return true;
}

double		
C3Matrix::Trace() const
{
	return m_a[0][0] + m_a[1][1] + m_a[2][2];
}

// (de)serialization
void
C3Matrix::Serialize( std::ostream& io_outstream ) const
{
	/*
	boost::archive::binary_oarchive oarch(io_outstream);
	oarch << *this;
	*/

	double tmp;
	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			// in_outstream << m_a[i][j] << " ";
			tmp = m_a[i][j];
			io_outstream.write((char*) &tmp, sizeof(tmp));
		}
	}
}

void
C3Matrix::Deserialize( std::istream& io_instream )
{
	/*
	boost::archive::binary_iarchive iarch(in_instream);
	iarch >> *this;
	*/

	double tmp;
	for(int j = 0; j < 3; j++)
	{
		for(int i = 0; i < 3; i++)
		{
			// in_instream >> m_a[i][j];
			io_instream.read((char*) &tmp, sizeof(tmp));
			m_a[i][j] = tmp;
		}
	}
}

// operators

C3Matrix 	
C3Matrix::operator+(const C3Matrix& in_mat) const
{
	C3Matrix result;
	for (int icol=0; icol < 3; icol++)
	{
		for (int jrow=0; jrow < 3; jrow++)
		{
			result.SetElement(icol, jrow, GetElement(icol, jrow) + in_mat.GetElement(icol, jrow));
		}
	}

	return result;
}

C3Matrix
C3Matrix::operator*(const C3Matrix& in_mat) const
{
	assert (m_a != NULL);

	C3Matrix result;
	C3Vector lhvec;
	C3Vector rhvec;
	for (int icol=0; icol < 3; icol++)
	{
		in_mat.GetColumn(icol, rhvec);
		for (int jrow=0; jrow < 3; jrow++)
		{
			GetRow(jrow, lhvec);
			double tmp = lhvec * rhvec;
			tmp = doubleEquals(tmp, 0.0, 0.000001) ? 0.0 : tmp;
			result.SetElement(icol, jrow, lhvec * rhvec);
		}
	}

	return result;
}

C3Vector
C3Matrix::operator*(const C3Vector& in_vec) const
{
	assert (m_a != NULL);

	C3Vector result;
	C3Vector lhvec;

	GetRow(0, lhvec);
	double v0 = lhvec * in_vec;
	v0 = doubleEquals(v0, 0.0, 0.000001) ? 0.0 : v0;

	GetRow(1, lhvec);
	double v1 = lhvec * in_vec;
	v1 = doubleEquals(v1, 0.0, 0.000001) ? 0.0 : v1;

	GetRow(2, lhvec);
	double v2 = lhvec * in_vec;
	v2 = doubleEquals(v2, 0.0, 0.000001) ? 0.0 : v2;

	result.Set(v0, v1, v2);

	return result;
}

void
C3Matrix::InitMatrixElements()
{
	m_a = new double* [3];
	for (int i=0; i < 3; i++)
		m_a[i] = new double [3];

	for (int i=0; i < 3; i++)
		for (int j=0; j < 3; j++)
			m_a[i][j] = 0.0;
}



















