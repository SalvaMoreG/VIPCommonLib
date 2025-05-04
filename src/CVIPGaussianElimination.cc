
#include "CVIPGaussianElimination.h"
#include "VIPconstants.h"

#include "CVIPMatrix.h"
#include "CVIPVector.h"


GaussianElimination::GaussianElimination( int in_dim )
	: m_ndim( in_dim )
	, m_matrixA(NULL)
	, m_b(NULL)
	, m_solutionP(NULL)
{
    m_matrixA = new CVIPMatrix( m_ndim, m_ndim );
    m_b = new CVIPVector( m_ndim, 0.0 );
	m_solutionP = new CVIPVector( m_ndim, 0.0 );
}

GaussianElimination::~GaussianElimination()
{
	delete m_matrixA;
	delete m_b;
	delete m_solutionP;
}

void
GaussianElimination::SetMatrixAElement( int irow, int icol, const double& in_value )
{
	m_matrixA->SetElement(irow, icol, in_value);
}

void
GaussianElimination::SetMatrixA( const CVIPMatrix& in_matrixA )
{
	(*m_matrixA) = in_matrixA;
}

void
GaussianElimination::SetVectorB( const CVIPVector& in_vectorB )
{
	(*m_b) = in_vectorB;
}

void
GaussianElimination::SolveSet( bool in_debug )
{
	// CREATE DIAGONAL FILE
	if (in_debug)
		cout << " =========== DIAGONAL FILE: ========== " << endl;
	for (int irow_sub = 0; irow_sub < m_ndim - 1; irow_sub++)
	{
	    int icol_diag = irow_sub;
		if (in_debug)
			cout << " =========== SUBTRACTING ROW " << irow_sub << endl;
        for (int krow = irow_sub+1; krow < m_ndim; krow++)
        {
            double factor(0.0);
            if (m_matrixA->GetElement(irow_sub, icol_diag) != 0)
                factor = m_matrixA->GetElement(krow, icol_diag)/m_matrixA->GetElement(irow_sub, icol_diag);
            SubtractRow(irow_sub, krow, factor);
        }
        if (in_debug)
			PrintEquations();
    }

    // NOW, THE OTHER WAY AROUND, SOLVE LAST ROW and MOVE UP.
	m_solutionP = new CVIPVector( m_ndim, 0.0 );
    if (m_matrixA->GetElement(m_ndim-1, m_ndim-1) == 0.0)
    {
        cout << "ERROR, no solution, m_matrixA.GetElement(m_ndim-1, m_ndim-1) = "
			 << m_matrixA->GetElement(m_ndim-1, m_ndim-1) << endl;
    }
    (*m_solutionP)[m_ndim-1] = (*m_b)[m_ndim-1]/m_matrixA->GetElement(m_ndim-1, m_ndim-1);

	if (in_debug)
		cout << " =========== SOLVING: ========== " << endl;
	for (int krow = m_ndim-2; krow >= 0; krow--)
	{
	    int kcol = krow;
        for (int irow_sub = m_ndim-1; irow_sub > krow; irow_sub--)
        {
			int icol_diag = irow_sub;
        	double factor(0.0);
        	if (m_matrixA->GetElement(irow_sub, icol_diag) != 0)
        		factor = m_matrixA->GetElement(krow, icol_diag)/m_matrixA->GetElement(irow_sub, icol_diag);
        	SubtractRow(irow_sub, krow, factor);
		}

        (*m_solutionP)[krow] = (*m_b)[krow]/m_matrixA->GetElement(krow, kcol);

		if (in_debug)
		{
			cout << " =========== SOLVING ROW: " << krow << endl;
			PrintEquations();
		}
    }

    if (in_debug)
	{
		cout << "--- FINAL SOLUTION:" << endl;
		for (int irow = 0; irow < m_ndim; irow++)
		{
			cout << (*m_solutionP)[irow] << endl;
		}
	}
}

void
GaussianElimination::GetSolutions( CVIPVector& out_solutionP )
{
	out_solutionP = (*m_solutionP);
}

void
GaussianElimination::PrintEquations()
{
	cout << " --- CURRENT SET OF EQUATIONS: " << endl;
	for (int irow = 0; irow < m_ndim; irow++)
	{
	    for (int icol = 0; icol < m_ndim; icol++)
		    cout << m_matrixA->GetElement( irow, icol) << " ";
		cout << " | " << (*m_b)[irow] << endl;
	}
	cout << endl;
}

void
GaussianElimination::SubtractRow(int irow_sub, int krow_from, const double& in_factor)
{
	(*m_b)[krow_from] -= in_factor * (*m_b)[irow_sub];
	for (int icol = 0; icol < m_ndim; icol++)
	{
	    double tmp = m_matrixA->GetElement( krow_from, icol) - in_factor * m_matrixA->GetElement( irow_sub, icol);
		m_matrixA->SetElement( krow_from, icol, tmp );
	}
}



