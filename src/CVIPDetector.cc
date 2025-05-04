
#include "CVIPDetector.h"

#include "CVIPRandom.h"
#include <cassert>

CDetector::CDetector()
	: m_dettype(DETTYPE_NONE)
	, m_eThreshold(0)
{
	m_size.Set(0, 0, 0);
	m_position.Set(0, 0, 0);;
	m_detVoxelSize.Set(0, 0, 0);
}

CDetector::~CDetector()
{
}

void
CDetector::SetType(DETECTOR_TYPE in_dettype)
{
	m_dettype = in_dettype;
}

CDetector::CDetector(const CDetector& in_obj)
{
	*this = in_obj;
}

CDetector&
CDetector::operator= (const CDetector& in_obj)
{
	m_size = in_obj.GetSize();
	m_position = in_obj.GetPosition();
    // return the existing object so we can chain this operator
    return *this;
}


bool
CDetector::IsInsideDetector( const C3Vector& in_position ) const
{
	double xmin = m_position.GetX() - 0.5 * m_size.GetX();
	double ymin = m_position.GetY() - 0.5 * m_size.GetY();
	double zmin = m_position.GetZ() - 0.5 * m_size.GetZ();

	/*
	bool isin = ( 	in_position.GetX() >= xmin && in_position.GetX() <= xmin + m_size.GetX()
		&& in_position.GetY() >= ymin && in_position.GetY() <= ymin + m_size.GetY()
		&& in_position.GetZ() >= zmin && in_position.GetZ() <= zmin + m_size.GetZ() );
	cout << "Isin: " << isin << " in_position: " << in_position
		 << "  lower: " << xmin << " " << ymin << " " << zmin
         << "  upper: " << xmin+m_size.GetX() << " " << ymin+m_size.GetY() << " " << zmin+m_size.GetZ() << endl;
	*/

	return ( 	in_position.GetX() >= xmin && in_position.GetX() <= xmin + m_size.GetX()
		&& in_position.GetY() >= ymin && in_position.GetY() <= ymin + m_size.GetY()
		&& in_position.GetZ() >= zmin && in_position.GetZ() <= zmin + m_size.GetZ() );
}

C3Vector
CDetector::CreateRandomRealHitPositionInDetector() const
{
	CVIPRandomUniform uniformDist(0.0, 1.0);
	double X = uniformDist.GetNewValue() * GetSize().GetX();
	double Y = uniformDist.GetNewValue() * GetSize().GetY();
	double Z = uniformDist.GetNewValue() * GetSize().GetZ();

	X = X - 0.5 * GetSize().GetX() + GetPosition().GetX();
	Y = Y - 0.5 * GetSize().GetY() + GetPosition().GetY();
	Z = Z - 0.5 * GetSize().GetZ() + GetPosition().GetZ();

	assert ( IsInsideDetector( C3Vector(X, Y, Z) ) );

	return C3Vector(X, Y, Z);
}

C3Vector
CDetector::GetVoxelizedHitPosition( const C3Vector& in_realHitPosition ) const
{
	// cout << "CreateDetPos in_realHitPos: " << in_realHitPosition << endl;
	assert ( IsInsideDetector( in_realHitPosition ) );

	double detVoxelSizeX = m_detVoxelSize.GetX();
	double detVoxelSizeY = m_detVoxelSize.GetY();
	double detVoxelSizeZ = m_detVoxelSize.GetZ();

	if (detVoxelSizeX == 0 && detVoxelSizeY == 0 && detVoxelSizeZ == 0) return in_realHitPosition;

	double x = in_realHitPosition.GetX();
	double y = in_realHitPosition.GetY();
	double z = in_realHitPosition.GetZ();

	// Make a grid of voxels, twice the size of the geometry resolution
	if (detVoxelSizeX > 0)
	{
		double xmin = - 0.5 * GetSize().GetX();		// size of entire detector
		x = x - xmin;
		assert (x >= 0);
		int BinNrX =  int( x / detVoxelSizeX );
		// cout << "bin nr: " << BinNrX << endl;
		x = xmin + ( BinNrX + 0.5 ) * detVoxelSizeX;
	}
	if (detVoxelSizeY > 0)
	{
		double ymin = - 0.5 * GetSize().GetY();
		y = y - ymin;
		assert (y >= 0);
		int BinNrY =  int( y / detVoxelSizeY );
		// cout << "bin nr: " << BinNrY << endl;
		y = ymin + ( BinNrY + 0.5 ) * detVoxelSizeY;
	}
	if (detVoxelSizeZ > 0)
	{
		double zmin = - 0.5 * GetSize().GetZ();
		z = z - zmin;
		assert (z >= 0);
		int BinNrZ =  int( z / detVoxelSizeZ );
		// cout << "bin nr: " << BinNrZ << endl;
		z = zmin + ( BinNrZ + 0.5 ) * detVoxelSizeZ;
	}
	// cout << "after applying voxelsize: " << x << " " << y << " " << z << endl;
	return C3Vector(x, y, z);
}









