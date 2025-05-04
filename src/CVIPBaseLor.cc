
#include "CVIPBaseLor.h"
#include "CVIPUtils.h"

#include <iostream>
#include <cassert>

using namespace std;

ostream& operator<<(ostream& os, const CVIPBaseLor& in_lor)
{
	os << "hit-1: " << in_lor.GetPosition1() << " E: " << in_lor.GetE1() << " ";
	os << "hit-2: " << in_lor.GetPosition2() << " E: " << in_lor.GetE2() << " ";
    return os;
}

CVIPBaseLor::CVIPBaseLor()
	: m_position1(0.0, 0.0, 0.0)
	, m_position2(0.0, 0.0, 0.0)
	, m_E1(0.0)
	, m_E2(0.0)
	, m_lorOrientation(LORDIR_NONE)
{
}

CVIPBaseLor::CVIPBaseLor(  const C3Vector& in_position1, const double& in_e1
                , const C3Vector& in_position2, const double& in_e2 )
	: m_position1(0.0, 0.0, 0.0)
	, m_position2(0.0, 0.0, 0.0)
	, m_lorOrientation(LORDIR_NONE)
{
    Set( in_position1, in_e1, in_position2, in_e2 );
}

CVIPBaseLor::~CVIPBaseLor()
{
}

int
CVIPBaseLor::Set(  const C3Vector& in_position1, const double& in_e1
        , const C3Vector& in_position2, const double& in_e2 )
{
    // GET MAX GRADIENT
	C3Vector LORVector = in_position2 - in_position1;

    double delX = fabs(LORVector.GetX());
    double delY = fabs(LORVector.GetY());
    double delZ = fabs(LORVector.GetZ());
	
	bool doswitch(false);
	if ( delZ > delX && delZ > delY )
	{
		m_lorOrientation = LORDIR_Z;		// z
		m_position1 = (in_position1.GetZ() < in_position2.GetZ()) ? in_position1 : in_position2;
		m_position2 = (in_position1.GetZ() < in_position2.GetZ()) ? in_position2 : in_position1;
		if (in_position1.GetZ() >= in_position2.GetZ()) doswitch = true;
	}
	else if ( delY > delX && delY > delZ )
	{
		m_lorOrientation = LORDIR_Y;		// y
		m_position1 = (in_position1.GetY() < in_position2.GetY()) ? in_position1 : in_position2;
		m_position2 = (in_position1.GetY() < in_position2.GetY()) ? in_position2 : in_position1;
		if (in_position1.GetY() >= in_position2.GetY()) doswitch = true;
	}
	else // if ( delX > delY && delX > delZ )
	{
		m_lorOrientation = LORDIR_X; 	// x
		m_position1 = (in_position1.GetX() < in_position2.GetX()) ? in_position1 : in_position2;
		m_position2 = (in_position1.GetX() < in_position2.GetX()) ? in_position2 : in_position1;
		if (in_position1.GetX() >= in_position2.GetX()) doswitch = true;
	}
	
	if (doswitch)
	{
		m_E1 = in_e2;
		m_E2 = in_e1;		
	}
	else
	{
		m_E1 = in_e1;
		m_E2 = in_e2;
	}
	
	return 0;
}

void
CVIPBaseLor::FindFovBins( const CVIPFieldOfView& in_fieldOfView
		, std::vector<unsigned int>& io_fovidxVec
		, std::vector<float>* io_fovweightsVec ) const
{
	double factor1pmil = GetMinFactorPromille( in_fieldOfView );
	/*
	cout << "factor1pmil: " << factor1pmil << endl;
	*/
	C3Vector lowerEdge( in_fieldOfView.GetLowerEdge() );
	C3Vector upperEdge( in_fieldOfView.GetUpperEdge() );
    
    C3Vector direction = GetPosition2() - GetPosition1();
    
    // <<<<<<<<<< Machiel, when hits lie inside FOV -- 2020-07-03
    C3Vector voxelSize = in_fieldOfView.GetVoxelSize();
    if (m_lorOrientation == LORDIR_X && direction.GetX() > 0)
    {
        if (GetPosition1().GetX() > lowerEdge.GetX())
            lowerEdge.SetX( GetPosition1().GetX() + voxelSize.GetX() );
        if (GetPosition2().GetX() < upperEdge.GetX())
            upperEdge.SetX( GetPosition2().GetX() - voxelSize.GetX() );
    }
    else if (m_lorOrientation == LORDIR_X)
    {
        if (GetPosition2().GetX() > lowerEdge.GetX())
            lowerEdge.SetX( GetPosition2().GetX() + voxelSize.GetX() );

        if (GetPosition1().GetX() < upperEdge.GetX())
            upperEdge.SetX( GetPosition1().GetX() - voxelSize.GetX() );
    }
    
    else if (m_lorOrientation == LORDIR_Y && direction.GetY() > 0 )
    {
        if (GetPosition1().GetY() > lowerEdge.GetY())
            lowerEdge.SetY( GetPosition1().GetY() + voxelSize.GetY() );
        if (GetPosition2().GetY() < upperEdge.GetY())
            upperEdge.SetY( GetPosition2().GetY() - voxelSize.GetY() );
    }
    else if (m_lorOrientation == LORDIR_Y)
    {
        if (GetPosition2().GetY() > lowerEdge.GetY())
            lowerEdge.SetY( GetPosition2().GetY() + voxelSize.GetY() );
        if (GetPosition1().GetY() < upperEdge.GetY())
            upperEdge.SetY( GetPosition1().GetY() - voxelSize.GetY() );
    }
    
    else if (m_lorOrientation == LORDIR_Z && direction.GetZ() > 0 )
    {
        if (GetPosition1().GetZ() > lowerEdge.GetZ())
            lowerEdge.SetZ( GetPosition1().GetZ() + voxelSize.GetZ() );
        if (GetPosition2().GetZ() < upperEdge.GetZ())
            upperEdge.SetZ( GetPosition2().GetZ() - voxelSize.GetZ() );
    }
    else if (m_lorOrientation == LORDIR_Z)
    {
        if (GetPosition2().GetZ() > lowerEdge.GetZ())
            lowerEdge.SetZ( GetPosition2().GetZ() + voxelSize.GetZ() );
        if (GetPosition1().GetZ() < upperEdge.GetZ())
            upperEdge.SetZ( GetPosition1().GetZ() - voxelSize.GetZ() );
    }
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
	double factor = 0.0;
	double factor_exit = 0.0;

	bool okay = RayBoxIntersection(GetPosition1(), direction, lowerEdge, upperEdge, factor, factor_exit);
	if ( !okay || factor <= 0 )
	{
		return;
	}

	// To be sure to not start outside of bounds (because of double == double comparisons...): 
	C3Vector currentPosition = GetPosition1() + direction * factor * (1 + factor1pmil);

	C3Vector lorInsideFov = direction * (factor_exit - factor);

	C3Vector minBounds, maxBounds;
	int index = -1;
	C3Vector prevPos;
	double weight;
	double sumweight = 0.0;
	//
	while (    currentPosition.GetX() <= upperEdge.GetX()
			&& currentPosition.GetY() <= upperEdge.GetY()
			&& currentPosition.GetZ() <= upperEdge.GetZ()
			&& in_fieldOfView.IsInBounds( currentPosition ) )
	{
		// if you have been here already, go a bit further....
		if ( in_fieldOfView.GetVoxelIndex( currentPosition ) == index )
		{
			// cout << "so, factor1pmil: " << factor1pmil << " so, factor: " << factor * (1 + factor1pmil) << " ";
			// cout << " old curPos: " << currentPosition << "  ";
			currentPosition = GetPosition1() + direction * factor * (1 + factor1pmil);

			// cout << "new curPos: " << currentPosition << " direction: " << direction 
			// 		<< " factor: " << factor << " f': " << factor * (1 + factor1pmil) << endl;
			// if ( in_fieldOfView.IsInBounds( currentPosition ) )
			// 	cout << " NEW new current index: " << in_fieldOfView.GetVoxelIndex( currentPosition ) << endl;
			// else 
			// 	cout << " NEW new current pos is out of bounds " << endl;

			if (    in_fieldOfView.IsInBounds( currentPosition )
				 && in_fieldOfView.GetVoxelIndex( currentPosition ) == index )
			{
				cout << "AGAIN!!! GONNA STOP LOOKING FOR MORE" << endl;
				return;
			}		
		}

		// deal with new position
		if ( in_fieldOfView.IsInBounds( currentPosition ) )
		{
			index = in_fieldOfView.GetVoxelIndex( currentPosition );

			// <<<<<<<<<<<<<<
			// int ix, iy, iz;
			// in_fieldOfView.GetVoxelIndices( index, ix, iy, iz);
			// cout << "factor: " << factor << "    pos:" << currentPosition << "   ";
			// cout << "with index: " << index << "   [" << ix << "," << iy << "," << iz << "]" << endl;
			// >>>>>>>>>>>>>>

			io_fovidxVec.push_back( (unsigned int) index );
	
			in_fieldOfView.GetVoxelBoundaries(index, minBounds, maxBounds);
			factor = NextRayIntersection( GetPosition1(), direction, minBounds, maxBounds );

			prevPos = currentPosition;
			currentPosition = GetPosition1() + direction * factor;

			if (io_fovweightsVec)
			{
				weight = (currentPosition - prevPos).GetLength()/direction.GetLength();
				// weight = (currentPosition - prevPos).GetLength()/lorInsideFov.GetLength();
				if ( weight < 0.0 || weight > 1.0 )
				{
					cout << "ERROR, weight: " << weight << endl;
					cout << "prevPos: " << prevPos << " currentPos: " << currentPosition 
						 << " LORdir: " << direction << endl;
					cout << " lorInsideFov: " << lorInsideFov << " with len: " << lorInsideFov.GetLength() << endl;
				}
				assert( weight >= 0.0 && weight <= 1.0 );

				io_fovweightsVec->push_back( (float) (weight) );

				sumweight += weight;
			}
		}

		// <<<<<<<<<<<<
		// cout << "io_fovidxVec size: " << io_fovidxVec.size() << endl;
		// cout << "origin: " << m_position1 << " direction: " << direction << " factor: " << factor << endl;
		// cout << index << " curPos: " << currentPosition << endl;
		// 
		// if (io_fovidxVec.size() > 500) exit(1);
		// >>>>>>>>>>>>>>>>>>>>
	}
}

/*
 GetMinFactorPromille()
 Input: 
	- CVIPFieldOfView object
 Output:
	- a small promilleage of the factor per voxel-size
 Method: 
   We have to move a promille (or less) just over the "edge" in order not to "fall back" into the previous bin..
*/

double 
CVIPBaseLor::GetMinFactorPromille( const CVIPFieldOfView& in_fieldOfView ) const
{
	C3Vector lowerEdge( in_fieldOfView.GetLowerEdge() );
	C3Vector upperEdge( in_fieldOfView.GetUpperEdge() );
	C3Vector direction = GetPosition2() - GetPosition1();

	double deltaFacXPerVoxel = (direction.GetX() > 0) ?
			  (upperEdge.GetX() - lowerEdge.GetX())/(direction.GetX() * in_fieldOfView.GetNumVoxelsX())
			: (lowerEdge.GetX() - upperEdge.GetX())/(direction.GetX() * in_fieldOfView.GetNumVoxelsX());
	double deltaFacYPerVoxel = (direction.GetY() > 0) ?
			  (upperEdge.GetY() - lowerEdge.GetY())/(direction.GetY() * in_fieldOfView.GetNumVoxelsY())
			: (lowerEdge.GetY() - upperEdge.GetY())/(direction.GetY() * in_fieldOfView.GetNumVoxelsY());	
	double deltaFacZPerVoxel = (direction.GetZ() > 0) ?
			  (upperEdge.GetZ() - lowerEdge.GetZ())/(direction.GetZ() * in_fieldOfView.GetNumVoxelsZ())
			: (lowerEdge.GetZ() - upperEdge.GetZ())/(direction.GetZ() * in_fieldOfView.GetNumVoxelsZ());
/*
cout << "deltaFacXPerVoxel: " << deltaFacXPerVoxel 
	 << " deltaFacYPerVoxel: " << deltaFacYPerVoxel
	 << " deltaFacZPerVoxel: " << deltaFacZPerVoxel << endl;
*/
	double factor1pmil = 0.0;

	if (fabs(deltaFacXPerVoxel) < fabs(deltaFacYPerVoxel) && fabs(deltaFacXPerVoxel) < fabs(deltaFacZPerVoxel) ) 
		factor1pmil = 0.001 * deltaFacXPerVoxel;
	else if (fabs(deltaFacYPerVoxel) <= fabs(deltaFacXPerVoxel) && fabs(deltaFacYPerVoxel) < fabs(deltaFacZPerVoxel) ) 
		factor1pmil = 0.001 * deltaFacYPerVoxel;
	// else if (fabs(deltaFacZPerVoxel) < fabs(deltaFacXPerVoxel) && fabs(deltaFacZPerVoxel) < fabs(deltaFacYPerVoxel) ) 
	else 
		factor1pmil = 0.001 * deltaFacZPerVoxel;

	return factor1pmil;
}

