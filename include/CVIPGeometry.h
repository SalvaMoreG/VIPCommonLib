
#ifndef _VIP_GEOMETRY_H__
#define _VIP_GEOMETRY_H__

#include "CVIP3Vector.h"
#include "CVIPFieldOfView.h"
#include "VIPconstants.h"
#include <sstream>

class CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipGeometryBase( const C3Vector& in_position );
	virtual				~CVipGeometryBase();

                    			CVipGeometryBase(const CVipGeometryBase& in_obj);
	virtual CVipGeometryBase& 	operator= (const CVipGeometryBase& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const {};

	virtual void 		SetPosition( const C3Vector in_position );
	const C3Vector&		GetPosition() const;

	virtual bool		IsInsideRegion( const C3Vector& in_position ) const = 0;
	virtual C3Vector 	GetMinimumPosition() const = 0;
	virtual C3Vector 	GetMaximumPosition() const = 0;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const = 0;

private:
	C3Vector			m_position;		// CENTER-POSITION OF OBJECT
};

// =======================================================================

class CVipBox : public CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipBox(  const C3Vector& in_position, const C3Vector& in_size);
	virtual				~CVipBox();

                    	CVipBox(const CVipBox& in_obj);
	virtual CVipBox& 	operator= (const CVipBox& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const
						{
							out_string << m_size;
						}

	bool				IsInsideRegion( const C3Vector& in_position ) const;
	virtual C3Vector 	GetMinimumPosition() const;
	virtual C3Vector 	GetMaximumPosition() const;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const;

private:

	void				SetSize( const C3Vector& in_size );
	const C3Vector&		GetSize() const;

	C3Vector			m_size;
};

// =======================================================================

class CVipSphere : public CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipSphere(  const C3Vector& in_position, const double& in_radius );
	virtual				~CVipSphere();
                   	
							CVipSphere(const CVipSphere& in_obj);
	virtual CVipSphere& 	operator= (const CVipSphere& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const
						{
							out_string << m_radius;
						}

	bool				IsInsideRegion( const C3Vector& in_position ) const;
	virtual C3Vector 	GetMinimumPosition() const;
	virtual C3Vector 	GetMaximumPosition() const;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const;

private:
	void				SetRadius( const double& in_radius );
	double				GetRadius() const;

	double 				m_radius;
};

/*
class CVipHemiSphere : public CVipSphere
{
public:
					// ctor and dtor
					CVipHemiSphere(  const C3Vector& in_position, const double& in_radius ) 
						: CVipSphere( in_position, in_radius ) {}
	virtual			~CVipHemiSphere() {};


};
*/

// =======================================================================

class CVipTube : public CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipTube(  const C3Vector& in_position, const double& in_radius, const double& in_height );
	virtual				~CVipTube();
                   	
						CVipTube(const CVipTube& in_obj);
	virtual CVipTube& 	operator= (const CVipTube& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const
						{
							out_string << m_radius << " " << m_height;
						}

	bool				IsInsideRegion( const C3Vector& in_position ) const;
	virtual C3Vector 	GetMinimumPosition() const;
	virtual C3Vector 	GetMaximumPosition() const;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const;

private:
	void				SetRadius( const double& in_radius );
	double				GetRadius() const;

	void				SetHeight( const double& in_height );
	double				GetHeight() const;

	double 				m_radius;
	double 				m_height;
};


// =======================================================================

class CVipCircle : public CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipCircle(  const C3Vector& in_position, const double& in_radius);
	virtual				~CVipCircle();

                    	CVipCircle(const CVipCircle& in_obj);
	virtual CVipCircle& operator= (const CVipCircle& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const
						{
							out_string << m_radius;
						}

	bool				IsInsideRegion( const C3Vector& in_position ) const;
	virtual C3Vector 	GetMinimumPosition() const;
	virtual C3Vector 	GetMaximumPosition() const;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const;

private:
	void				SetRadius( const double& in_radius );
	double				GetRadius() const;

	double 				m_radius;
};

// =======================================================================

// WARNING: Lines are always horizontal (along X).... TODO: make it more general.
class CVipLine : public CVipGeometryBase
{
public: 
						// ctor and dtor
						CVipLine(  const C3Vector& in_position, const double& in_length);
	virtual				~CVipLine();

                    	CVipLine(const CVipLine& in_obj);
	virtual CVipLine&   operator= (const CVipLine& in_obj);


						// get and set
	virtual double		GetVolume() const;
	virtual void		PrintProperty(std::stringstream& out_string) const
						{
							out_string << m_length;
						}

	bool				IsInsideRegion( const C3Vector& in_position ) const;
	virtual C3Vector 	GetMinimumPosition() const;
	virtual C3Vector 	GetMaximumPosition() const;
	virtual bool		IsVoxelInsideGeometry( int in_voxelindex, const CVIPFieldOfView& in_fov ) const;

private:
	void				SetLength( const double& in_length );
	double				GetLength() const;

	double 				m_length;
};


#endif

