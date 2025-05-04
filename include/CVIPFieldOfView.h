
#ifndef _C_FIELDOFVIEW__H__
#define _C_FIELDOFVIEW__H__

#include "CVIP3Vector.h"
#include "CVIP3Matrix.h"

// useful class
class CVIPFieldOfView
{
public: 
	// constructor and destructor
							CVIPFieldOfView();
	virtual					~CVIPFieldOfView();

	// copy and assign constructors
							CVIPFieldOfView(const CVIPFieldOfView& in_obj);
	CVIPFieldOfView& 		operator= (const CVIPFieldOfView& in_obj);


	void					Initialize( const std::string& io_geometryFileName );

    // Set
    void					SetFieldOfView(   int in_bins_x, const double& in_xmin, const double& in_xmax
											, int in_bins_y, const double& in_ymin, const double& in_ymax
											, int in_bins_z, const double& in_zmin, const double& in_zmax);

	// Get	
	int						GetNumberOfVoxels() const;
	inline int				GetNumVoxelsX() const { return m_nVoxelsX; }
	inline int				GetNumVoxelsY() const { return m_nVoxelsY; }
	inline int				GetNumVoxelsZ() const { return m_nVoxelsZ; }

	inline const C3Vector&	GetLowerEdge() const { return m_edgeMin; }
	inline const C3Vector&	GetUpperEdge() const { return m_edgeMax; }

	inline C3Vector		 	GetSize() const 	  	{ return (m_edgeMax - m_edgeMin); }
	inline const C3Vector& 	GetVoxelSize() const 	{ return m_voxelSize; }
	double 					GetVoxelVolume() const;

	bool 					IsValidIndex(int in_index) const;
	bool					AreValidIndices(int in_ix, int in_iy, int in_iz) const;
	
	bool					GetVoxelCentre(int in_index, C3Vector& io_coords) const;
	bool					GetVoxelCentre(int in_ix, int in_iy, int in_iz
										, double& io_x, double& io_y, double& io_z) const;
	bool					GetVoxelBoundaries(int in_index, C3Vector& io_min, C3Vector& io_max) const;
	bool					IsInVoxelBounds( int in_index, const C3Vector& in_coords ) const;

	int						GetVoxelIndex(const C3Vector& in_coords ) const;    // returns -1 if not inside image
	bool					GetVoxelIndices(const C3Vector& in_coords, int& io_ix, int& io_iy, int& io_iz ) const;

	int						GetVoxelIndex( int in_ix, int in_iy, int in_iz ) const;
	bool					GetVoxelIndices(int in_index, int& io_ix, int& io_iy, int& io_iz) const;

	bool					IsInBounds( const C3Vector& coords ) const;
	C3Vector				GetCenterCoord() const;
	C3Vector				GetVoxelCornerCoordinates(int in_fovIndex, int cornerIndex) const;

	// WARNING!!! The next function should be used with extreme care!!!
	// NOTE!!! Only works with angles of +- 90 degrees and +- 180 degrees!!! (i.e. cos and sin = +- 1 and 0)
	void					Rotate( const C3Matrix& in_matrix );
	
	bool					IsSymmetric() const;

private:
	void					ReadSetupParameters( const std::string& io_geometryFileName );

private:
    
	// define output stream operator
	friend std::ostream& 		operator<<( std::ostream& os, const CVIPFieldOfView& fov);
    
	// edges
	C3Vector 				m_edgeMin;
	C3Vector 				m_edgeMax;

	// number of voxels in each direction
	int 					m_nVoxelsX;
	int 					m_nVoxelsY;
	int 					m_nVoxelsZ;

	// size
	C3Vector				m_voxelSize;
};

#endif
