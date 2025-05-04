
#ifndef _VIP_GAUSSIAN_ELIMINATION_H__
#define _VIP_GAUSSIAN_ELIMINATION_H__

#include "VIPconstants.h"

#include "CVIPMatrix.h"
#include "CVIPVector.h"

class GaussianElimination
{
public:
				GaussianElimination( int in_dim );
				~GaussianElimination();

	void		SetMatrixAElement( int irow, int icol, const double& in_value );
	void		SetMatrixA( const CVIPMatrix& in_matrixA );
	void		SetVectorB( const CVIPVector& in_vectorB );

	void		SolveSet( bool in_debug );
	void		GetSolutions( CVIPVector& out_solutionP );

	void 		PrintEquations();

private:

	// Subtracting Rows
	void 		SubtractRow(int irow_sub, int krow_from, const double& in_factor);

	// Data
	int 		m_ndim;

	CVIPMatrix* m_matrixA;
	CVIPVector* m_b;
	CVIPVector* m_solutionP;
};


#endif
