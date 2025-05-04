
#include "CVIPGeometry.h"

#include <algorithm>

#include <iostream>
#include <cassert>

using namespace std;

/*
CVipGeometryBase::CVipGeometryBase() 
	: m_position(0, 0, 0) 
{
}
*/

CVipGeometryBase::CVipGeometryBase( const C3Vector& in_position ) 
	: m_position(in_position) 
{
}

CVipGeometryBase::~CVipGeometryBase() 
{
}

CVipGeometryBase::CVipGeometryBase(const CVipGeometryBase& in_obj)
{
	m_position = in_obj.GetPosition();
	// cout << "Copy CTOR BASE, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
}

CVipGeometryBase& 	
CVipGeometryBase::operator= (const CVipGeometryBase& in_obj)
{
	m_position = in_obj.GetPosition();
	// cout << "Assign op BASE, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
	return *this;
}


double	
CVipGeometryBase::GetVolume() const
{
	cout << "NOT DEFINED" << endl;
	return 0.0;
}

void
CVipGeometryBase::SetPosition( const C3Vector in_position ) 
{ 
	m_position = in_position; 
}

const C3Vector&	
CVipGeometryBase::GetPosition() const 
{ 
	return m_position; 
}

// ===========================================================================

CVipBox::CVipBox(  const C3Vector& in_position, const C3Vector& in_size) 
	: CVipGeometryBase( in_position )
{ 
	SetSize(in_size); 
}

CVipBox::~CVipBox()
{
}

CVipBox::CVipBox(const CVipBox& in_obj)
	: CVipGeometryBase(in_obj)
{
	m_size = in_obj.GetSize();
	// cout << "Copy CTOR CVipBox, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
}
	
CVipBox&
CVipBox::operator= (const CVipBox& in_obj)
{
	CVipGeometryBase::operator=(in_obj);
	m_size = in_obj.GetSize();
	// cout << "Assign op CVipBox, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
	return *this;
}

double
CVipBox::GetVolume() const
{
	// cout << "box volume, using size: " << m_size << endl;
	return m_size.GetX() * m_size.GetY() * m_size.GetZ();
}

bool
CVipBox::IsInsideRegion( const C3Vector& in_position ) const
{
	C3Vector minPos = GetMinimumPosition();
	C3Vector maxPos = GetMaximumPosition();
	if (       in_position.GetX() > minPos.GetX() 
			&& in_position.GetX() < maxPos.GetX() 
			&& in_position.GetY() > minPos.GetY() 
			&& in_position.GetY() < maxPos.GetY()
			&& in_position.GetZ() > minPos.GetZ() 
			&& in_position.GetZ() < maxPos.GetZ() )
	{
		return true;
	} 
	return false;
}

bool
CVipBox::IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const
{
	C3Vector voxel_edgeMin;
	C3Vector voxel_edgeMax;
	in_fov.GetVoxelBoundaries( in_voxelindex, voxel_edgeMin, voxel_edgeMax);
			
	// Box: 
	bool IsInsideX = (   fabs(voxel_edgeMin.GetX() - GetPosition().GetX()) < 0.5 * m_size.GetX() 
					  || fabs(GetPosition().GetX() - voxel_edgeMax.GetX()) < 0.5 * m_size.GetX());
	bool IsInsideY = (   fabs(voxel_edgeMin.GetY() - GetPosition().GetY()) < 0.5 * m_size.GetY() 
					  || fabs(GetPosition().GetY() - voxel_edgeMax.GetY()) < 0.5 * m_size.GetY());
	bool IsInsideZ = (   fabs(voxel_edgeMin.GetZ() - GetPosition().GetZ()) < 0.5 * m_size.GetZ() 
					  || fabs(GetPosition().GetZ() - voxel_edgeMax.GetZ()) < 0.5 * m_size.GetZ());

	return IsInsideX && IsInsideY && IsInsideZ;
}

C3Vector 	
CVipBox::GetMinimumPosition() const
{
	C3Vector aVector(  GetPosition() - (m_size*0.5) );
	return aVector;
}

C3Vector
CVipBox::GetMaximumPosition() const
{
	C3Vector aVector(  GetPosition() + (m_size*0.5) );
	return aVector;
}

const C3Vector&	
CVipBox::GetSize() const 
{ 
	return m_size; 
}

void
CVipBox::SetSize( const C3Vector& in_size )
{
	m_size = in_size;
}

// =====================================================================================

CVipSphere::CVipSphere(  const C3Vector& in_position, const double& in_radius ) 
	: CVipGeometryBase( in_position )
{ 
	SetRadius(in_radius); 
}

CVipSphere::~CVipSphere() 
{
}

CVipSphere::CVipSphere(const CVipSphere& in_obj)
	: CVipGeometryBase(in_obj)
{
	m_radius = in_obj.GetRadius();
	// cout << "Copy CTOR CVipSphere, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
}

CVipSphere& 	
CVipSphere::operator= (const CVipSphere& in_obj)
{
	SetPosition( in_obj.GetPosition() );
	m_radius = in_obj.GetRadius();
	// cout << "Assign op CVipSphere, in_obj volume: " << in_obj.GetVolume() << " this vol: " << GetVolume() << endl;
	return *this;
}


double
CVipSphere::GetVolume() const
{
	// cout << "sphere volume, using radius: " << m_radius << endl;
	return (4.0/3.0) * kPI * m_radius*m_radius*m_radius;
}

bool
CVipSphere::IsInsideRegion( const C3Vector& in_position ) const
{
	C3Vector relativePos (in_position - GetPosition());
	double distance = relativePos.GetLength();
	if ( distance < m_radius )
	{
		return true;
	} 
	return false;
}

bool
CVipSphere::IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const
{
	C3Vector voxel_edgeMin;
	C3Vector voxel_edgeMax;
	in_fov.GetVoxelBoundaries( in_voxelindex, voxel_edgeMin, voxel_edgeMax);

// cout << "[" << in_voxelindex << "] voxel_edgeMin: " << voxel_edgeMin << " max: " << voxel_edgeMax ;
	//	Sphere:
	/*
	double x1 = voxel_edgeMin.GetX() - GetPosition().GetX();
	double x2 = voxel_edgeMax.GetX() - GetPosition().GetX();
	double xso = std::min(x1, x2); 
	cout << "e.g. X: min: " << x1 << " max: " << x2 << " so: " << xso << endl;
	double minDistX = std::min( fabs(voxel_edgeMin.GetX() - GetPosition().GetX())
							  , fabs(voxel_edgeMax.GetX() - GetPosition().GetX()) );
	double minDistY = std::min( fabs(voxel_edgeMin.GetY() - GetPosition().GetY())
 							  , fabs(voxel_edgeMax.GetY() - GetPosition().GetY()) );
	double minDistZ = std::min( fabs(voxel_edgeMin.GetZ() - GetPosition().GetZ())
							  , fabs(voxel_edgeMax.GetZ() - GetPosition().GetZ()) );
	*/

	double minDistX(0.0), minDistY(0.0), minDistZ(0.0); 

	if (voxel_edgeMin.GetX() > GetPosition().GetX())
		minDistX = voxel_edgeMin.GetX() - GetPosition().GetX();
	else if (voxel_edgeMax.GetX() < GetPosition().GetX())
		minDistX = GetPosition().GetX() - voxel_edgeMax.GetX();

	if (voxel_edgeMin.GetY() > GetPosition().GetY())
		minDistY = voxel_edgeMin.GetY() - GetPosition().GetY();
	else if (voxel_edgeMax.GetY() < GetPosition().GetY())
		minDistY = GetPosition().GetY() - voxel_edgeMax.GetY();

	if (voxel_edgeMin.GetZ() > GetPosition().GetZ())
		minDistZ = voxel_edgeMin.GetZ() - GetPosition().GetZ();
	else if (voxel_edgeMax.GetZ() < GetPosition().GetZ())
		minDistZ = GetPosition().GetZ() - voxel_edgeMax.GetZ();

	assert (minDistX >= 0);
	assert (minDistY >= 0);
	assert (minDistZ >= 0);

	double dist = sqrt(minDistX*minDistX + minDistY*minDistY + minDistZ*minDistZ);

// cout << " dist_vector: [" << minDistX << "," << minDistY << "," << minDistZ << "] dist: " << dist << endl;

	return (dist < m_radius);
}

C3Vector 	
CVipSphere::GetMinimumPosition() const
{
	C3Vector aVector(  GetPosition() - (C3Vector(1, 1, 1) * m_radius) );
	return aVector;
}

C3Vector
CVipSphere::GetMaximumPosition() const
{
	C3Vector aVector(  GetPosition() + (C3Vector(1, 1, 1) * m_radius) );
	return aVector;
}

double
CVipSphere::GetRadius() const 
{ 
	return m_radius; 
}

void
CVipSphere::SetRadius( const double& in_radius )
{
	m_radius = in_radius;
}

// =====================================================================================

CVipTube::CVipTube(  const C3Vector& in_position, const double& in_radius, const double& in_height ) 
	: CVipGeometryBase( in_position )
{ 
	SetRadius(in_radius); 
	SetHeight(in_height); 
}

CVipTube::~CVipTube() 
{
}

CVipTube::CVipTube(const CVipTube& in_obj)
	: CVipGeometryBase(in_obj)
{
	m_radius = in_obj.GetRadius();
	m_height = in_obj.GetHeight();
}

CVipTube& 	
CVipTube::operator= (const CVipTube& in_obj)
{
	SetPosition( in_obj.GetPosition() );
	m_radius = in_obj.GetRadius();
	m_height = in_obj.GetHeight();
	return *this;
}


double
CVipTube::GetVolume() const
{
	return kPI * m_radius*m_radius * m_height;
}

bool
CVipTube::IsInsideRegion( const C3Vector& in_position ) const
{
	C3Vector relativePos (in_position - GetPosition());
	double distanceXY = sqrt(relativePos.GetX()*relativePos.GetX() + relativePos.GetY()*relativePos.GetY());
	double distanceZ = relativePos.GetZ();
	if ( fabs(distanceZ) < 0.5 * m_height && distanceXY < m_radius )
	{
		return true;
	} 
	return false;
}

bool
CVipTube::IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const
{
return false;
	// TODO: figure this out
/*
	C3Vector voxel_edgeMin;
	C3Vector voxel_edgeMax;
	in_fov.GetVoxelBoundaries( in_voxelindex, voxel_edgeMin, voxel_edgeMax);

	C3Vector voxelCentre;
	in_fov.GetVoxelCentre(in_voxelindex, voxelCentre);
	double minDistX(0.0), minDistY(0.0), minDistZ(0.0); 

	if (voxel_edgeMin.GetX() > GetPosition().GetX())
		minDistX = voxel_edgeMin.GetX() - GetPosition().GetX();
	else if (voxel_edgeMax.GetX() < GetPosition().GetX())
		minDistX = GetPosition().GetX() - voxel_edgeMax.GetX();

	if (voxel_edgeMin.GetY() > GetPosition().GetY())
		minDistY = voxel_edgeMin.GetY() - GetPosition().GetY();
	else if (voxel_edgeMax.GetY() < GetPosition().GetY())
		minDistY = GetPosition().GetY() - voxel_edgeMax.GetY();

	if (voxel_edgeMin.GetZ() > GetPosition().GetZ())
		minDistZ = voxel_edgeMin.GetZ() - GetPosition().GetZ();
	else if (voxel_edgeMax.GetZ() < GetPosition().GetZ())
		minDistZ = GetPosition().GetZ() - voxel_edgeMax.GetZ();

	assert (minDistX >= 0);
	assert (minDistY >= 0);
	assert (minDistZ >= 0);

	double dist = sqrt(minDistX*minDistX + minDistY*minDistY + minDistZ*minDistZ);

	return (dist < m_radius);
*/
}

C3Vector 	
CVipTube::GetMinimumPosition() const
{
	C3Vector aVector(  GetPosition() - C3Vector(m_radius, m_radius, m_height) );
	return aVector;
}

C3Vector
CVipTube::GetMaximumPosition() const
{
	C3Vector aVector(  GetPosition() + C3Vector(m_radius, m_radius, m_height) );
	return aVector;
}

double
CVipTube::GetRadius() const 
{ 
	return m_radius; 
}

void
CVipTube::SetRadius( const double& in_radius )
{
	m_radius = in_radius;
}

double
CVipTube::GetHeight() const 
{ 
	return m_height; 
}

void
CVipTube::SetHeight( const double& in_height )
{
	m_height = in_height;
}

// ============================================================================================

CVipCircle::CVipCircle( const C3Vector& in_position, const double& in_radius ) 
	: CVipGeometryBase( in_position )
{ 
	SetRadius(in_radius); 
}

CVipCircle::~CVipCircle() 
{
}

CVipCircle::CVipCircle(const CVipCircle& in_obj)
	: CVipGeometryBase(in_obj)
{
	m_radius = in_obj.GetRadius();
}

CVipCircle& 	
CVipCircle::operator= (const CVipCircle& in_obj)
{
	SetPosition( in_obj.GetPosition() );
	m_radius = in_obj.GetRadius();

	return *this;
}


double
CVipCircle::GetVolume() const
{
	return 0.0;
}

bool
CVipCircle::IsInsideRegion( const C3Vector& in_position ) const
{
	C3Vector relativePos (in_position - GetPosition());
	relativePos.SetZ(0.0);
	double distance = relativePos.GetLength();
	if (   distance < m_radius 
        && doubleEquals(in_position.GetZ(), GetPosition().GetZ(), 0.0001) )
	{
/*
cout << "in_position.GetZ(): " << in_position.GetZ() << " " 
     << " GetPosition().GetZ(): " << GetPosition().GetZ() << " " 
	 << " equal?: " << doubleEquals(in_position.GetZ(), GetPosition().GetZ(), 0.0001)  << endl;
*/
		return true;
	} 
	return false;
}

bool
CVipCircle::IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const
{
	C3Vector voxel_edgeMin;
	C3Vector voxel_edgeMax;
	in_fov.GetVoxelBoundaries( in_voxelindex, voxel_edgeMin, voxel_edgeMax);

	C3Vector voxelCentre;
	in_fov.GetVoxelCentre(in_voxelindex, voxelCentre);
	double minDistX(0.0), minDistY(0.0); 

	if (voxel_edgeMin.GetX() > GetPosition().GetX())
		minDistX = voxel_edgeMin.GetX() - GetPosition().GetX();
	else if (voxel_edgeMax.GetX() < GetPosition().GetX())
		minDistX = GetPosition().GetX() - voxel_edgeMax.GetX();

	if (voxel_edgeMin.GetY() > GetPosition().GetY())
		minDistY = voxel_edgeMin.GetY() - GetPosition().GetY();
	else if (voxel_edgeMax.GetY() < GetPosition().GetY())
		minDistY = GetPosition().GetY() - voxel_edgeMax.GetY();

	assert (minDistX >= 0);
	assert (minDistY >= 0);

	double dist = sqrt(minDistX*minDistX + minDistY*minDistY);

	// return (    dist < m_radius 
    //          && doubleEquals(voxelCentre.GetZ(), GetPosition().GetZ(), 0.0001) );

	double voxSizeZ = in_fov.GetVoxelSize().GetZ();
	bool isIn = (   dist < m_radius
			&& fabs(voxelCentre.GetZ() - GetPosition().GetZ()) < voxSizeZ );

	return isIn;
}

C3Vector 	
CVipCircle::GetMinimumPosition() const
{
	C3Vector aVector(  GetPosition() - (C3Vector(1, 1, 0) * m_radius) );
	return aVector;
}

C3Vector
CVipCircle::GetMaximumPosition() const
{
	C3Vector aVector(  GetPosition() + (C3Vector(1, 1, 0) * m_radius) );
	return aVector;
}

double
CVipCircle::GetRadius() const 
{ 
	return m_radius; 
}

void
CVipCircle::SetRadius( const double& in_radius )
{
	m_radius = in_radius;
}

// =====================================================================================

CVipLine::CVipLine( const C3Vector& in_position, const double& in_length ) 
	: CVipGeometryBase( in_position )
{ 
	SetLength(in_length); 
}

CVipLine::~CVipLine() 
{
}

CVipLine::CVipLine(const CVipLine& in_obj)
	: CVipGeometryBase(in_obj)
{
	m_length = in_obj.GetLength();
}

CVipLine& 	
CVipLine::operator= (const CVipLine& in_obj)
{
	SetPosition( in_obj.GetPosition() );
	m_length = in_obj.GetLength();
	return *this;
}


double
CVipLine::GetVolume() const
{
	return 0.0;
}

bool
CVipLine::IsInsideRegion( const C3Vector& in_position ) const
{
	double distance = fabs( in_position.GetX() - GetPosition().GetX() );
	if (    distance < 0.5 * m_length 
		 && doubleEquals(in_position.GetY(), GetPosition().GetY(), 0.0001)
		 && doubleEquals(in_position.GetZ(), GetPosition().GetZ(), 0.0001) )
	{
		return true;
	} 
	return false;
}

bool
CVipLine::IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const
{
	C3Vector voxel_edgeMin;
	C3Vector voxel_edgeMax;
	in_fov.GetVoxelBoundaries( in_voxelindex, voxel_edgeMin, voxel_edgeMax);

	double minDistX(0.0); 
	if (voxel_edgeMin.GetX() > GetPosition().GetX())
		minDistX = voxel_edgeMin.GetX() - GetPosition().GetX();
	else if (voxel_edgeMax.GetX() < GetPosition().GetX())
		minDistX = GetPosition().GetX() - voxel_edgeMax.GetX();

	assert (minDistX >= 0);

	C3Vector voxelCentre;
	in_fov.GetVoxelCentre(in_voxelindex, voxelCentre);

	/*
	return (   minDistX < 0.5 * m_length 
			&& doubleEquals(voxelCentre.GetY(), GetPosition().GetY(), 0.0001)
			&& doubleEquals(voxelCentre.GetZ(), GetPosition().GetZ(), 0.0001) );
	*/

	double voxSizeY = in_fov.GetVoxelSize().GetY();
	double voxSizeZ = in_fov.GetVoxelSize().GetZ();

	bool isIn = (   minDistX < 0.5 * m_length 
			&& fabs(voxelCentre.GetY() - GetPosition().GetY()) < voxSizeY
			&& fabs(voxelCentre.GetZ() - GetPosition().GetZ()) < voxSizeZ );
	// cout << "CVipLine::IsVoxelInsideGeometry: " << fck << endl;
	// cout << "Line @:" << GetPosition() << " min: " << GetMinimumPosition() << " max: " << GetMaximumPosition() << endl;

	// C3Vector coords;
	// in_fov.GetVoxelCentre(in_voxelindex, coords);
	// cout << "        " << in_voxelindex << " @: " << coords << endl;

	return isIn;
}

C3Vector 	
CVipLine::GetMinimumPosition() const
{
//	C3Vector delta = (C3Vector(1, 0, 0) * 0.5 * m_length);
	C3Vector aVector( GetPosition() - (C3Vector(1, 0, 0) * 0.5 * m_length) );
// cout << "MinPos, pos: " << GetPosition() << " delta: " << delta << " minPos: " << aVector << " len: " << m_length << endl;
	return aVector;
}

C3Vector
CVipLine::GetMaximumPosition() const
{
	C3Vector aVector(  GetPosition() + (C3Vector(1, 0, 0) * 0.5 * m_length) );
	return aVector;
}

double
CVipLine::GetLength() const 
{ 
	return m_length; 
}

void
CVipLine::SetLength( const double& in_length )
{
	m_length = in_length;
}

