
#pragma once

#ifndef _VIPWILDERMAN_UTILS_H__
#define _VIPWILDERMAN_UTILS_H__

#include "CVIP3Vector.h"
#include "CVIPUtils.h"
#include "CVIPFieldOfView.h"

#include <vector>
#include <cassert>
#include <cstdlib>

enum SEARCH_DIRECTION
{
	SEARCH_NONE = 0,  
	SEARCH_NEXT_X,  
	SEARCH_PREV_X,
	SEARCH_NEXT_Y,  
	SEARCH_PREV_Y
};

// ------------------ struct EllipsParameters ---------------
struct EllipsParameters
{
	double A, B, C, D, E, F;
	double axis_a, axis_b, centre_x, centre_y;
};


// ------------------ class Comparator ---------------------
class Comparator 
{
public: 
	enum COMPARISON 
	{
		COMPARE_NONE = 0, 
		COMPARE_SMALLER,
		COMPARE_BIGGER
	};
	
					Comparator( COMPARISON in_comparison ) : m_comparison( in_comparison ) {}
					~Comparator() {}

	inline bool 	Compare(const double& in_var1, const double& in_var2) const
					{
						if (m_comparison == COMPARE_SMALLER)
							return in_var1 < in_var2;
						else if (m_comparison == COMPARE_BIGGER)
							return in_var1 > in_var2;
						else
						{
							assert( false );
							return false;
						}
					}

private:
	// copy and assign constructors not defined, cannot be used
					Comparator(const Comparator& in_obj);
	Comparator& 	operator= (const Comparator& in_obj);
	
	// data
	COMPARISON 		m_comparison;
};

// ------------------ various methods -----------------------
void
PrintSearchDirection( SEARCH_DIRECTION in_searchDirection );

int WildermanMarch( const double& in_slicePositionZ   
			, const C3Vector& in_coneOrigin, const C3Vector& in_coneAxis, const double& in_coneAngle
			, const CVIPFieldOfView& in_fieldOfView, std::vector<unsigned int>& io_fovidxVec, bool doFillBins = true );

void 
March( const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
			, const C3Vector& in_intersect1, const C3Vector& in_intersect2, const CVIPFieldOfView& in_fieldOfView
			, std::vector<unsigned int>& io_fovidxVec, SEARCH_DIRECTION in_searchDirection );

void 
SubMarch( int in_fovidx, const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
			, const C3Vector& in_intersect1, const C3Vector& in_intersect2
			, const CVIPFieldOfView& in_fieldOfView, std::vector<unsigned int>& io_fovidxVec
			, SEARCH_DIRECTION in_searchDirection, unsigned int in_end_fovidx );

int
FindNearestPoint( SEARCH_DIRECTION in_searchDirection, const double& in_currentPos, 
		int in_numSolutions1, const C3Vector& in_sol1, const C3Vector& in_sol2, 
		int in_numSolutions2, const C3Vector& in_sol3, const C3Vector& in_sol4);

void
CheckEquals(int& io_numSolutions, const C3Vector& in_currentPos, C3Vector& io_sol1, C3Vector& io_sol2);

int 
CompareThreePoints( const double& in_currentPos, const double& in_sol1Pos, const double& in_sol2Pos
					, const Comparator& in_comparator );

void
FillBins( AXIS in_axis, int in_currentidxX, int in_currentidxY, int in_currentidxZ, int in_nextidxQ
         , const CVIPFieldOfView& in_fieldOfView, std::vector<unsigned int>& io_fovidxVec
         , bool& io_doStop, unsigned int in_end_fovidx );

int 
FindIntersectionsInFOV(const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
        , const C3Vector& in_Edge, AXIS in_toSolve
		, const C3Vector& in_minBounds, const C3Vector& in_maxBounds
		, C3Vector& io_intersect1, C3Vector& io_intersect2
		, QuadraticFactors& io_quadraticFactors );

void 
FromConeEquationToQuadraticFactor( const C3Vector& in_coneAxis, const C3Vector& in_scatPos, const double& in_coneAngle 
	, AXIS in_toSolve, const C3Vector& in_edge, QuadraticFactors& io_quadraticFactors );

/*
void
FromConeEquationToQuadraticFactorOld( const C3Vector& in_coneAxis, const C3Vector& in_scatPos, const double& in_coneAngle 
	, AXIS in_toSolve, const C3Vector& in_edge, QuadraticFactors& io_quadraticFactors );
*/
	
void
ShowFovVector( const CVIPFieldOfView& in_fieldOfView, const std::vector<unsigned int>& in_fovidxVec );

void
CalculateEllipsFactors( const double& in_slicePositionZ   
			, const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_comptonAngle
			, const CVIPFieldOfView& in_fieldOfView, EllipsParameters& io_ellipsParameters );


void
AdditionalParameters( EllipsParameters& io_ellipsParameters );

#endif

