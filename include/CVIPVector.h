#ifndef _CVIPCOMMONLIB_Vector__H__
#define _CVIPCOMMONLIB_Vector__H__

#pragma once

#include <iostream>
#include <fstream>
#include <ostream>

using namespace std;

class CVIPVector
{
public: 
	// constructor and destructor
					CVIPVector( int nvoxels );
					CVIPVector( int nvoxels, const double& in_initvalue );
	virtual			~CVIPVector();

	// copy and assign constructors
					CVIPVector(const CVIPVector& in_obj);
	CVIPVector& 	operator= (const CVIPVector& in_obj);

	// Get/Set methods
	int				GetSize() const 		{ return m_nvoxels; }

	// operators
	double&			operator[](int index);
	const double&	operator[](int index) const;

	CVIPVector 		operator*(const double& in_coeff) const;
	double 			operator*(const CVIPVector&) const;    

	CVIPVector 		operator+(const CVIPVector&) const;
	CVIPVector 		operator-(const CVIPVector&) const;
   
	bool			operator==(const CVIPVector& in_obj) const;
	bool			operator!=(const CVIPVector& in_obj) const;

private:
	// define output stream operator
	friend ostream& 	operator<<(ostream& os, const CVIPVector& vec);

	double*			m_data;
	int				m_nvoxels;
};


#endif
