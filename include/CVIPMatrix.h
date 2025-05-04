#ifndef _VIPCOMMONLIB_MATRIX_H___
#define _VIPCOMMONLIB_MATRIX_H___

#include "CVIPVector.h"

#include <iostream>
#include <fstream>
#include <ostream>

class CVIPMatrix
{
public: 
	// constructor and destructor
					CVIPMatrix( int in_nrows, int in_ncolumns );
                    CVIPMatrix( int in_nrows, int in_ncolumns, const double& init_value );
	virtual			~CVIPMatrix();

	// copy and assign constructors
					CVIPMatrix(const CVIPMatrix& in_obj);
	CVIPMatrix& 	operator= (const CVIPMatrix& in_obj);

	//
	int 			GetNumberOfRows() const { return m_nrows; }
	int 			GetNumberOfColumns() const { return m_ncolumns; }

	// 
	const double&	GetElement( int irow, int icolumn ) const;
	void			SetElement( int irow, int icolumn, const double& in_value );

	// operators
	CVIPVector		operator*(const CVIPVector&) const;
    CVIPMatrix      operator+(const CVIPMatrix& in_m) const;
    CVIPMatrix      operator-(const CVIPMatrix& in_m) const;
   
	bool			operator==(const CVIPMatrix& in_obj) const;
	bool			operator!=(const CVIPMatrix& in_obj) const;

private:
	// aux methods
	void			InitMatrixElements( const double& init_value = 0.0 );

	// define output stream operator
	friend ostream& 	operator<<(ostream& os, const CVIPMatrix& vec);

	double**		m_data;
	int				m_nrows;
	int				m_ncolumns;

};

#endif
