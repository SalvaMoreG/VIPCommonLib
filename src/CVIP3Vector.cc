
#include "CVIP3Vector.h"
#include "CVIP3Matrix.h"
#include "CVIPUtils.h"

#include <cmath>

const double kEPSILON = 0.000001;

using namespace std;

ostream& operator<<(ostream& os, const C3Vector& vec)
{
    os << vec.GetX() << ' ' << vec.GetY() << ' ' << vec.GetZ();
    return os;
}

istream& operator>>( istream &is, C3Vector& vec )
{ 
    is >> vec.m_x >> vec.m_y >> vec.m_z;
    return is;
}

C3Vector::C3Vector()
	: m_x(0.0)
	, m_y(0.0)
	, m_z(0.0)
{
}

C3Vector::C3Vector(const double& in_x, const double& in_y, const double& in_z)
	: m_x(in_x)
	, m_y(in_y)
	, m_z(in_z)
{
}

C3Vector::~C3Vector()
{
}

// copy constructor
C3Vector::C3Vector(const C3Vector& in_obj)
{
	*this = in_obj;
}

// assignment operator
C3Vector&
C3Vector::operator= (const C3Vector& in_obj)
{
	if (this != &in_obj)
	{
		m_x = in_obj.GetX();
		m_y = in_obj.GetY();
		m_z = in_obj.GetZ();
	}
	return *this;
}

bool
C3Vector::operator==(const C3Vector &in_obj) const
{
	return (   doubleEquals( m_x, in_obj.GetX(), kEPSILON )
			&& doubleEquals( m_y, in_obj.GetY(), kEPSILON )
			&& doubleEquals( m_z, in_obj.GetZ(), kEPSILON ) );
}

bool
C3Vector::operator!=(const C3Vector &in_obj) const
{
    return !(*this == in_obj);
}

void
C3Vector::Set (const C3Vector& in_obj)
{
	*this = in_obj;
}

void
C3Vector::Set (const double& in_x, const double& in_y, const double& in_z)
{
	m_x = in_x;
	m_y = in_y;
	m_z = in_z;
}

double
C3Vector::GetLength() const
{
	double len = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
	return len;
}

double
C3Vector::Get2DProjectionLength( AXIS in_perpendicularAxis ) const
{
	if ( in_perpendicularAxis  == Axis_Z )
		return sqrt(m_x*m_x + m_y*m_y);
	else if ( in_perpendicularAxis  == Axis_Y )
		return sqrt(m_x*m_x + m_z*m_z);
	else if ( in_perpendicularAxis  == Axis_X )
		return sqrt(m_z*m_z + m_y*m_y);
	else
		return 0.0;
}

void
C3Vector::Normalize(double in_newlen)
{
	double len = GetLength();
	if (len > 0)
	{
		(*this) = (*this) * (in_newlen / len);
	}
}

double
C3Vector::GetScalarProductCosAngle(const C3Vector& in_vec) const
{
    double cosw = (*this) * in_vec;
	double len1 = GetLength();
	double len2 = in_vec.GetLength();

	cosw = ((len1 * len2) != 0) ? (cosw / (len1 * len2)) : 0;
    return cosw;
}

double
C3Vector::GetScalarProductAngleRadians(const C3Vector& in_vec) const
{
	double cosangle = GetScalarProductCosAngle(in_vec);
    double test(0);
    if ( fabs(cosangle) < 1)
        test = acos(cosangle);

    //  cout << "cosangle: " << cosangle << " so test(rad): " << test << endl;
	return test;
}

// (De)Serialization
void
C3Vector::Serialize( std::ostream& io_outstream ) const
{
	io_outstream.write((char*) &m_x, sizeof(m_x));
	io_outstream.write((char*) &m_y, sizeof(m_y));
	io_outstream.write((char*) &m_z, sizeof(m_z));
}

void
C3Vector::Deserialize( std::istream& io_instream )
{
	io_instream.read((char*) &m_x, sizeof(m_x));
	io_instream.read((char*) &m_y, sizeof(m_y));
	io_instream.read((char*) &m_z, sizeof(m_z));
}

// operators

// + (plus) operator
C3Vector
C3Vector::operator+ (const C3Vector& in_v) const
{
      C3Vector result;
      result.m_x = (m_x + in_v.GetX());
      result.m_y = (m_y + in_v.GetY());
      result.m_z = (m_z + in_v.GetZ());
      return result;
}

// - (minus) operator
C3Vector
C3Vector::operator- (const C3Vector& in_v) const
{
      C3Vector result;
      result.m_x = (m_x - in_v.GetX());
      result.m_y = (m_y - in_v.GetY());
      result.m_z = (m_z - in_v.GetZ());
      return result;
}

// inproduct
double
C3Vector::operator*(const C3Vector& in_v) const
{
	double result = m_x * in_v.GetX() + m_y * in_v.GetY() + m_z * in_v.GetZ();
	return result;
}

C3Vector
C3Vector::operator*(const C3Matrix& in_matrix) const
{
	C3Vector result;
	C3Vector colvec;
	in_matrix.GetColumn(0, colvec);
	result.SetX( *this * colvec);
	in_matrix.GetColumn(1, colvec);
	result.SetY( *this * colvec);
	in_matrix.GetColumn(2, colvec);
	result.SetZ( *this * colvec);

	return result;
}

C3Vector
C3Vector::operator*(double const& in_coeff) const
{
	C3Vector result;
	result.m_x = (m_x * in_coeff);
	result.m_y = (m_y * in_coeff);
	result.m_z = (m_z * in_coeff);
	return result;
}

// outer-product
C3Matrix
C3Vector::OuterProduct(const C3Vector& in_v) const
{
    C3Matrix matrix;
    // void		SetElement(int col, int row, double value);
    matrix.SetElement(0, 0, m_x * in_v.GetX());
    matrix.SetElement(1, 0, m_x * in_v.GetY());
    matrix.SetElement(2, 0, m_x * in_v.GetZ());

    matrix.SetElement(0, 1, m_y * in_v.GetX());
    matrix.SetElement(1, 1, m_y * in_v.GetY());
    matrix.SetElement(2, 1, m_y * in_v.GetZ());

    matrix.SetElement(0, 2, m_z * in_v.GetX());
    matrix.SetElement(1, 2, m_z * in_v.GetY());
    matrix.SetElement(2, 2, m_z * in_v.GetZ());

	return matrix;
}














