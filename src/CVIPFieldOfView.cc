
#include "CVIPFieldOfView.h"
#include "VIPconstants.h"
#include "CVIPERRORCODES.h"

#include <cassert>
#include <cstdlib>

using namespace std;

// =================================================================================

ostream& operator<<(ostream& os, const CVIPFieldOfView& fov)
{
    os << fov.GetNumVoxelsX() << ' ' << fov.GetNumVoxelsY() << ' ' << fov.GetNumVoxelsZ() << ' ' 
       << fov.GetLowerEdge() << ' ' << fov.GetUpperEdge();
    return os;    
}

// =================================================================================

CVIPFieldOfView::CVIPFieldOfView()
 : m_edgeMin(0, 0, 0)
 , m_edgeMax(0, 0, 0)
 , m_nVoxelsX(0)
 , m_nVoxelsY(0)
 , m_nVoxelsZ(0)
{
}

CVIPFieldOfView::~CVIPFieldOfView()
{
}

// copy and assign constructors
CVIPFieldOfView::CVIPFieldOfView(const CVIPFieldOfView& in_obj)
{
	*this = in_obj;
}

CVIPFieldOfView&
CVIPFieldOfView::operator= (const CVIPFieldOfView& in_obj)
{
	m_nVoxelsX = in_obj.GetNumVoxelsX();
	m_nVoxelsY = in_obj.GetNumVoxelsY();
	m_nVoxelsZ = in_obj.GetNumVoxelsZ();

	m_edgeMin = in_obj.GetLowerEdge();
	m_edgeMax = in_obj.GetUpperEdge();

	m_voxelSize = in_obj.GetVoxelSize();

	/*
	C3Vector test( (m_edgeMax.GetX() - m_edgeMin.GetX())/m_nVoxelsX
			  	 , (m_edgeMax.GetY() - m_edgeMin.GetY())/m_nVoxelsY
				 , (m_edgeMax.GetZ() - m_edgeMin.GetZ())/m_nVoxelsZ ); 
	assert( test == m_voxelSize );
	*/

	return *this;
}

void
CVIPFieldOfView::Initialize( const std::string& io_geometryFileName )
{
	ReadSetupParameters( io_geometryFileName );
}

void
CVIPFieldOfView::SetFieldOfView(   int in_bins_x, const double& in_xmin, const double& in_xmax
								 , int in_bins_y, const double& in_ymin, const double& in_ymax
								 , int in_bins_z, const double& in_zmin, const double& in_zmax )
{
	assert( in_bins_x > 0 && in_bins_y > 0 && in_bins_z > 0 );
	assert( in_xmax > in_xmin && in_ymax > in_ymin && in_zmax > in_zmin );

	m_nVoxelsX = in_bins_x;
	m_nVoxelsY = in_bins_y;
	m_nVoxelsZ = in_bins_z;

	m_edgeMin.Set( in_xmin, in_ymin, in_zmin );
	m_edgeMax.Set( in_xmax, in_ymax, in_zmax );

	m_voxelSize.Set( (in_xmax - in_xmin)/in_bins_x
				   , (in_ymax - in_ymin)/in_bins_y
				   , (in_zmax - in_zmin)/in_bins_z ); 
}

void
CVIPFieldOfView::ReadSetupParameters( const std::string& io_geometryFileName )
{
	ifstream conffile(io_geometryFileName.c_str(), ios::in);
	if (!conffile.is_open())
	{
		cout << "WARNING!!! Configure file not open: " << io_geometryFileName << endl;
		throw VIPERRORCODES::VIP_ERROR_FOV_FILE_NOT_FOUND;
	}

	std::string tmpStr, tmpDum;
	double tmpVal;

	while (!conffile.eof())
    {
		conffile >> tmpStr >> tmpVal >> tmpDum;
		if (conffile.fail() && !conffile.eof())
		{
			cout << "ERROR! FORMAT ERROR in parameters file" << endl;
			exit(1);
		}
		if (!conffile.eof())
		{
            if (tmpStr == "minX")
			{
				m_edgeMin.SetX(tmpVal);
			}
			else if (tmpStr == "minY")
			{
				m_edgeMin.SetY(tmpVal);
			}
			else if (tmpStr == "minZ")
			{
				m_edgeMin.SetZ(tmpVal);
			}
            else if (tmpStr == "maxX")
			{
				m_edgeMax.SetX(tmpVal);
			}
			else if (tmpStr == "maxY")
			{
				m_edgeMax.SetY(tmpVal);
			}
			else if (tmpStr == "maxZ")
			{
				m_edgeMax.SetZ(tmpVal);
			}
			else if (tmpStr == "nVoxelsX")
			{
				m_nVoxelsX = (int) (tmpVal);
			}
			else if (tmpStr == "nVoxelsY")
			{
				m_nVoxelsY = (int) (tmpVal);
			}
			else if (tmpStr == "nVoxelsZ")
			{
				m_nVoxelsZ = (int) (tmpVal);
			}
			else
			{
				cout << "ERROR. Unknown parameter: " << tmpStr << endl;
				exit(1);
			}
		}
	}
	// some checks
	assert ( m_edgeMax.GetX() > m_edgeMin.GetX() );
	assert ( m_edgeMax.GetY() > m_edgeMin.GetY() );
	assert ( m_edgeMax.GetZ() > m_edgeMin.GetZ() );

	// establish voxel size
	m_voxelSize.Set( (m_edgeMax.GetX() - m_edgeMin.GetX())/m_nVoxelsX
			  	   , (m_edgeMax.GetY() - m_edgeMin.GetY())/m_nVoxelsY
				   , (m_edgeMax.GetZ() - m_edgeMin.GetZ())/m_nVoxelsZ ); 
}

// *******************************************************************************************************
int
CVIPFieldOfView::GetNumberOfVoxels() const
{
	return m_nVoxelsX * m_nVoxelsY * m_nVoxelsZ;
}

// *******************************************************************************************************

double
CVIPFieldOfView::GetVoxelVolume() const
{
	return GetVoxelSize().GetX() * GetVoxelSize().GetY() * GetVoxelSize().GetZ();
}

// *******************************************************************************************************

bool 
CVIPFieldOfView::IsValidIndex(int in_index) const
{
	return ( in_index >= 0 && in_index < GetNumberOfVoxels() );
}

bool
CVIPFieldOfView::AreValidIndices(int in_ix, int in_iy, int in_iz) const
{
	return (    in_ix >= 0 && in_ix < m_nVoxelsX 
	         && in_iy >= 0 && in_iy < m_nVoxelsY 
	         && in_iz >= 0 && in_iz < m_nVoxelsZ );
}

bool
CVIPFieldOfView::GetVoxelCentre(int in_index, C3Vector& io_coords) const
{
	if ( !IsValidIndex( in_index ) )
	{
		return false;
	}
	
	int ix = 0, iy = 0, iz = 0;
	if ( !GetVoxelIndices(in_index, ix, iy, iz) )
	{
		return false;
	}

	double x, y, z;
	GetVoxelCentre( ix, iy, iz, x, y, z);

    io_coords.Set(x, y, z);
	
	return true;
}

bool
CVIPFieldOfView::GetVoxelCentre(int in_ix, int in_iy, int in_iz
							  , double& io_x, double& io_y, double& io_z) const
{
	if ( !AreValidIndices(in_ix, in_iy, in_iz)  )
	{
		return false;
	}
	
    io_x = m_edgeMin.GetX() + (in_ix + 0.5) * m_voxelSize.GetX();
    io_y = m_edgeMin.GetY() + (in_iy + 0.5) * m_voxelSize.GetY();
    io_z = m_edgeMin.GetZ() + (in_iz + 0.5) * m_voxelSize.GetZ();
	return true;
}

bool
CVIPFieldOfView::GetVoxelBoundaries(int in_index, C3Vector& io_min, C3Vector& io_max) const
{
	if ( !IsValidIndex( in_index ) )
	{
		return false;
	}
	
	int ix = 0, iy = 0, iz = 0;
	GetVoxelIndices(in_index, ix, iy, iz);

	io_min.Set (  m_edgeMin.GetX() + ix * m_voxelSize.GetX()
				, m_edgeMin.GetY() + iy * m_voxelSize.GetY()
				, m_edgeMin.GetZ() + iz * m_voxelSize.GetZ());

	io_max = io_min + m_voxelSize;
	return true;
}

bool
CVIPFieldOfView::IsInVoxelBounds( int in_index, const C3Vector& in_coords ) const
{
	if ( !IsValidIndex( in_index ) )
	{
		return false;
	}
	
	C3Vector edgeMin, edgeMax;
	GetVoxelBoundaries( in_index, edgeMin, edgeMax);

	bool result =      (edgeMin.GetX() <= in_coords.GetX() && in_coords.GetX() <= edgeMax.GetX());
	result = result && (edgeMin.GetY() <= in_coords.GetY() && in_coords.GetY() <= edgeMax.GetY());
	result = result && (edgeMin.GetZ() <= in_coords.GetZ() && in_coords.GetZ() <= edgeMax.GetZ());

	return result;
}


C3Vector
CVIPFieldOfView::GetVoxelCornerCoordinates(int in_index, int in_cornerIndex) const
{
	assert( in_cornerIndex >= 0 && in_cornerIndex < 8 );
	int ix = 0, iy = 0, iz = 0;
	GetVoxelIndices(in_index, ix, iy, iz);

	// cornerIndex = 0, x = voxel_xmin, y = voxel_ymin, z = voxel_zmin
	// cornerIndex = 1, x = voxel_xmax, y = voxel_ymin, z = voxel_zmin
	// cornerIndex = 2, x = voxel_xmin, y = voxel_ymax, z = voxel_zmin
	// cornerIndex = 3, x = voxel_xmax, y = voxel_ymax, z = voxel_zmin
	// cornerIndex = 4, x = voxel_xmin, y = voxel_ymin, z = voxel_zmax
	// cornerIndex = 5, x = voxel_xmax, y = voxel_ymin, z = voxel_zmax
	// cornerIndex = 6, x = voxel_xmin, y = voxel_ymax, z = voxel_zmax
	// cornerIndex = 7, x = voxel_xmax, y = voxel_ymax, z = voxel_zmax

	double x = (in_cornerIndex % 2 == 0) ? m_edgeMin.GetX() + ix*m_voxelSize.GetX() 
									     : m_edgeMin.GetX() + (ix+1)*m_voxelSize.GetX();
	double y = (in_cornerIndex == 0 || in_cornerIndex == 1 || in_cornerIndex == 4 || in_cornerIndex == 5) 
						? m_edgeMin.GetY() + iy*m_voxelSize.GetY() : m_edgeMin.GetY() + (iy+1)*m_voxelSize.GetY();
	double z = (in_cornerIndex <= 3) ? m_edgeMin.GetZ() + iz*m_voxelSize.GetZ() : m_edgeMin.GetZ() + (iz+1)*m_voxelSize.GetZ();

	C3Vector coordinates(x, y, z);
	return coordinates;
}

// *******************************************************************************************************

int
CVIPFieldOfView::GetVoxelIndex(const C3Vector& in_coords ) const
{
	int ix = 0, iy = 0, iz = 0;

	if (!GetVoxelIndices(in_coords, ix, iy, iz ))
	{
		return -1;
	}

    int index = GetVoxelIndex( ix, iy, iz );

	assert( index >= 0 && index < GetNumberOfVoxels() );
	return index;
}

bool
CVIPFieldOfView::GetVoxelIndices(const C3Vector& in_coords, int& io_ix, int& io_iy, int& io_iz ) const
{
    assert ( IsInBounds( in_coords ) );
	if ( !IsInBounds( in_coords ) )
	{
		return false;
	}

    io_ix = (int) ((in_coords.GetX() - m_edgeMin.GetX()) / m_voxelSize.GetX());
    io_iy = (int) ((in_coords.GetY() - m_edgeMin.GetY()) / m_voxelSize.GetY());
    io_iz = (int) ((in_coords.GetZ() - m_edgeMin.GetZ()) / m_voxelSize.GetZ());

	// When on the upper edges, count as the previous index
	if (io_ix == m_nVoxelsX) io_ix = m_nVoxelsX - 1; 
	if (io_iy == m_nVoxelsY) io_iy = m_nVoxelsY - 1; 
	if (io_iz == m_nVoxelsZ) io_iz = m_nVoxelsZ - 1; 

	return true;
}


int
CVIPFieldOfView::GetVoxelIndex( int in_ix, int in_iy, int in_iz ) const
{
	if ( !AreValidIndices(in_ix, in_iy, in_iz)  )
	{
		return -1;
	}	
	return in_iz*m_nVoxelsX*m_nVoxelsY + in_iy*m_nVoxelsX + in_ix;
}


bool
CVIPFieldOfView::GetVoxelIndices(int in_index, int& io_ix, int& io_iy, int& io_iz) const
{
	if ( !IsValidIndex( in_index ) )
	{
		return false;
	}	
	
    io_iz = in_index / (m_nVoxelsX * m_nVoxelsY);
    int rest = in_index % (m_nVoxelsX * m_nVoxelsY);

    io_iy = rest / m_nVoxelsX;
    io_ix = rest % m_nVoxelsX;

	return true;
}


// *******************************************************************************************************

bool
CVIPFieldOfView::IsInBounds( const C3Vector& in_coords ) const
{
	bool result =      (m_edgeMin.GetX() <= in_coords.GetX() && in_coords.GetX() <= m_edgeMax.GetX());
	result = result && (m_edgeMin.GetY() <= in_coords.GetY() && in_coords.GetY() <= m_edgeMax.GetY());
	result = result && (m_edgeMin.GetZ() <= in_coords.GetZ() && in_coords.GetZ() <= m_edgeMax.GetZ());

	// <<<<<<<<<<<<<<<
	/*
	cout << "edgeMinX: " << m_edgeMin.GetX() << " incoordsX: " << in_coords.GetX() << endl;
	cout << "IsInB, a:  " << (m_edgeMin.GetX() <= in_coords.GetX()) << "  "; 
	cout << "IsInB, b:  " << (in_coords.GetX() <= m_edgeMax.GetX()) << "  "; 
	cout << "IsInB, c:  " << (m_edgeMin.GetY() <= in_coords.GetY()) << "  "; 
	cout << "IsInB, d:  " << (in_coords.GetY() <= m_edgeMax.GetY()) << "  "; 
	cout << "IsInB, e:  " << (m_edgeMin.GetZ() <= in_coords.GetZ()) << "  "; 
	cout << "IsInB, f:  " << (in_coords.GetZ() <= m_edgeMax.GetZ()) << endl;
	*/
	// >>>>>>>>>>>>>>>

	return result;
}

// *******************************************************************************************************

C3Vector
CVIPFieldOfView::GetCenterCoord() const
{		
	return (m_edgeMin + m_edgeMax)*0.5;
}

// *******************************************************************************************************

// NOTE!!! Only works with angles of +- 90 degrees and +- 180 degrees!!! (i.e. cos and sin = +- 1 and 0)
void
CVIPFieldOfView::Rotate( const C3Matrix& in_matrix )
{
	// 1. Edges
	m_edgeMin = in_matrix * m_edgeMin;
	m_edgeMax = in_matrix * m_edgeMax;

	// min values might be bigger than max values!!!!!!
	C3Vector tmpmin( std::min(m_edgeMin.GetX(), m_edgeMax.GetX())
				   , std::min(m_edgeMin.GetY(), m_edgeMax.GetY())
				   , std::min(m_edgeMin.GetZ(), m_edgeMax.GetZ()) );

	m_edgeMax.Set( std::max(m_edgeMin.GetX(), m_edgeMax.GetX())
				 , std::max(m_edgeMin.GetY(), m_edgeMax.GetY())
				 , std::max(m_edgeMin.GetZ(), m_edgeMax.GetZ()) );

	m_edgeMin = tmpmin;

	// cout << "AFTER, m_edgeMin: " << m_edgeMin << ", m_edgeMax: " << m_edgeMax << endl;

	// 2. Voxel Size
	C3Vector row;
	in_matrix.GetRow(0, row);
	double xsize = std::abs( row.GetX() * m_voxelSize.GetX() + row.GetY() * m_voxelSize.GetY() + row.GetZ() * m_voxelSize.GetZ() );

	in_matrix.GetRow(1, row);
	double ysize = std::abs( row.GetX() * m_voxelSize.GetX() + row.GetY() * m_voxelSize.GetY() + row.GetZ() * m_voxelSize.GetZ() );

	in_matrix.GetRow(2, row);
	double zsize = std::abs( row.GetX() * m_voxelSize.GetX() + row.GetY() * m_voxelSize.GetY() + row.GetZ() * m_voxelSize.GetZ() );

	m_voxelSize.Set( xsize, ysize, zsize );

	// cout << "AFTER, voxelSize: " << m_voxelSize << endl;

	// 3. Number of voxels
	m_nVoxelsX = (int) (  (m_edgeMax.GetX() - m_edgeMin.GetX())/m_voxelSize.GetX() + 0.5 );
	m_nVoxelsY = (int) (  (m_edgeMax.GetY() - m_edgeMin.GetY())/m_voxelSize.GetY() + 0.5 );
	m_nVoxelsZ = (int) (  (m_edgeMax.GetZ() - m_edgeMin.GetZ())/m_voxelSize.GetZ() + 0.5 );

	// cout << "AFTER, m_nVoxelsX: " << m_nVoxelsX << " m_nVoxelsY: " << m_nVoxelsY << " m_nVoxelsZ: " << m_nVoxelsZ << endl;
}

// *******************************************************************************************************

bool
CVIPFieldOfView::IsSymmetric() const
{
	return (  m_nVoxelsX == m_nVoxelsY && m_nVoxelsX == m_nVoxelsZ
		&& m_edgeMin.GetX() == -1.0 * m_edgeMax.GetX()
		&& m_edgeMin.GetY() == -1.0 * m_edgeMax.GetY()
		&& m_edgeMin.GetZ() == -1.0 * m_edgeMax.GetZ() );
}	
	

