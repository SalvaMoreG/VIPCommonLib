
#pragma once

#ifndef _VIPCOMMON_C3VECTOR___
#define _VIPCOMMON_C3VECTOR___

#include "CVIPUtils.h"

#include <iostream>
#include <fstream>

class C3Matrix;

class C3Vector
{
public:
				C3Vector();
				C3Vector(const double& in_x, const double& in_y, const double& in_z);
	virtual		~C3Vector();

	// copy and assign constructors
				C3Vector(const C3Vector& in_obj);
	C3Vector& 	operator= (const C3Vector& in_obj);

	// set function
	void			Set (const double& in_x, const double& in_y, const double& in_z);
	void			Set (const C3Vector& in_obj);
	inline double 	Get ( int index ) const
					{
						if ( index < 0 && index >= 3 )
						{
							std::cout << "ERROR, C3Vector::Get( index ), index: " << index 
									  << ", return 0.0" << std::endl;
							return 0.0;
						}
						if (index == 0) return m_x;
						else if (index == 1) return m_y;
						else if (index == 2) return m_z;

						return 0.0;
					}
	inline double	GetX() const { return m_x; }
	inline double	GetY() const { return m_y; }
	inline double	GetZ() const { return m_z; }
	inline void		SetX(const double& in_val) { m_x = in_val; }
	inline void		SetY(const double& in_val) { m_y = in_val; }
	inline void		SetZ(const double& in_val) { m_z = in_val; }

	// useful functions
	double		GetLength() const;
	double		Get2DProjectionLength( AXIS in_perpendicularAxis ) const;
	void		Normalize(double in_newlen = 1.0);
	double		GetScalarProductCosAngle(const C3Vector& in_vec) const;       	// scalarproduct 
		// NOTE!!!! Already divided by the lengths of both vectors!!!!!!
	double		GetScalarProductAngleRadians(const C3Vector& in_vec) const;       		// scalarproduct angle in radians

	// (De)Serialization
	void		Serialize( std::ostream& outstream ) const;
	void		Deserialize( std::istream& instream );

	// operators
	C3Vector 	operator+(const C3Vector&) const;       		// operator+()
	C3Vector 	operator-(const C3Vector&) const;       		// operator-()

	double 		operator*(const C3Vector& in_vector) const;     // 
	C3Vector 	operator*(const C3Matrix& in_matrix) const;     // 
	C3Vector 	operator*(const double& in_coeff) const;        // multiplication factor

	bool		operator==(const C3Vector& in_obj) const; 		// ==
	bool		operator!=(const C3Vector& in_obj) const;		// !=

	C3Matrix	OuterProduct(const C3Vector& in_v) const; 		// outer product

private:
	// data
	double 						m_x, m_y, m_z;

	// define output stream operator
	friend std::ostream& 		operator<<( std::ostream& os, const C3Vector& vec);
    friend std::istream&        operator>>( std::istream& is, C3Vector& vec );

	// BOOST SERIALIZATION STUFF
	/*
	friend class boost::serialization::access;
	// template functions should always be defined where they are defined (in header file)
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & m_x & m_y & m_z;
    }
	*/
};

#endif




