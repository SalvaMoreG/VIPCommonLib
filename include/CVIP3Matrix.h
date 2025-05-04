#pragma once

#ifndef _VIPCOMMON_C3MATRIX___
#define _VIPCOMMON_C3MATRIX___

#include "CVIP3Vector.h"
#include "CVIPUtils.h"

#include "CVIPUtils.h"

#include <iostream>
#include <fstream>

class C3Matrix
{
public:

	// ctors and dtor
				C3Matrix();
				C3Matrix(const C3Vector& in_col1, const C3Vector& in_col2, const C3Vector& in_col3);
				C3Matrix( AXIS in_axis, const double& in_angle);		// angles in radians

	virtual		~C3Matrix();


	// copy and assign constructors
				C3Matrix(const C3Matrix& in_obj);
	C3Matrix& 	operator= (const C3Matrix& in_obj);

	// set function
	void		Set (const C3Vector& in_col1, const C3Vector& in_col2, const C3Vector& in_col3);
	void		Set( AXIS in_axis, const double& in_angle);
	double		GetElement(int col, int row) const;
	void		SetElement(int col, int row, const double& value);
	void		SetUnity( const double& in_factor  );

	void 		GetColumn(int col, C3Vector& io_vec) const;
	void 		GetRow(int row, C3Vector& io_vec) const;

	bool		IsSymmetric() const;
	double		Trace() const;

	// useful functions
	// double		GetInverse() const;		// not implemented yet

	// serialization
	void		Serialize( std::ostream& outstream ) const;
	void		Deserialize( std::istream& instream );

	// operators
	C3Matrix 	operator+(const C3Matrix&) const;

	C3Matrix 	operator*(const C3Matrix&) const;
	C3Vector 	operator*(const C3Vector&) const;    
   
	bool		operator==(const C3Matrix& in_obj) const; 		// ==
	bool		operator!=(const C3Matrix& in_obj) const;		// !=

private:

	// data
	double** 				m_a;

	// initialize data
	void 					InitMatrixElements();

	// define output stream operator
	friend std::ostream& 	operator<<(std::ostream& os, const C3Matrix& vec);
};

#endif
