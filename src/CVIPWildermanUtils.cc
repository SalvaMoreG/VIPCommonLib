
#include "CVIPWildermanUtils.h"
#include "VIPconstants.h"

#include <cassert>
#include <cstdlib>
#include <algorithm>

using namespace std;

//	const bool m_debugOutput( true );
const bool m_debugOutput( false );

void
PrintSearchDirection( SEARCH_DIRECTION in_searchDirection )
{
	if (in_searchDirection == SEARCH_NONE)
		cout << " NONE ";
	else if (in_searchDirection == SEARCH_NEXT_X)
		cout << " NEXT X ";
	else if (in_searchDirection == SEARCH_PREV_X)
		cout << " PREV X ";
	else if (in_searchDirection == SEARCH_NEXT_Y)  
		cout << " NEXT Y ";
	else if (in_searchDirection == SEARCH_PREV_Y)
		cout << " PREV Y ";
	else
		assert( false );
}

// ****************************************************************************************************

// Based on the Wilderman article "Fast Algorithm for List Mode Back-Projection of Compton Scatter Camera Data"
//  ( EEE TRANSACTIONS ON NUCLEAR SCIENCE, VOL. 45, NO. 3 )
// In the Wilderman article, the scenario numbering scheme is based on starting with one intersection, how many intersections are there
// in a) the same edge b) the opposite edge) c) if nothing else, in the orthogonal edge
// My scenario numbering is more or less based on how many total intersections there are,
//		before sub-dividing it on the configuration of these intersections

/*
#intersections               Wilderman paper scenario             My code scenario
0 (cone outside FOV)         -                                    0
0 (cone inside FOV)          1                                    1

1,3,5,7                      "never happens"                      ignored

2 (orthogonal edges)         2                                    2
2 (opposite edges)           3                                    2
2 (only one edge)            5                                    5

4 (2 in orthogal edges)      5                                    3
4 (1 orthogonal/1 opposite)  4                                    4
4 (2 in opposite edge)       6                                    6
4 (1 in each edge)           3                                    11

4 (2 in same orthogonal edge) -                                   10 (very very unlikely or maybe impossible?)

6 (4 in orthogonal edges)    one of the above                     7
8 (2 in each edge)           one of the above                     8
else (???)                   one of the above (???)               9

*/


int WildermanMarch( const double& in_slicePositionZ   
			, const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_comptonAngle
			, const CVIPFieldOfView& in_fieldOfView, std::vector<unsigned int>& io_fovidxVec, bool doFillBins )
{
	int scenario = -1;
	io_fovidxVec.clear();

	// TODO:
	// Note: this is only for the case that the detector is oriented in the XY-plane
	//		 and the source is away from the detector along the Z-axis
	//
	// 1. Find Solutions for All Edges
	const C3Vector& minBounds = in_fieldOfView.GetLowerEdge();
	const C3Vector& maxBounds = in_fieldOfView.GetUpperEdge();

	double dummy = 0.0;
	int numSolutionsEdge[4];
	C3Vector edge_sol1[4], edge_sol2[4];

	QuadraticFactors quadraticFactors;

	C3Vector edge( minBounds.GetX(), dummy, in_slicePositionZ );	
	if (m_debugOutput) cout << "Y_axis, left-edge: " << edge << " "; 
	numSolutionsEdge[0] = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, edge, Axis_Y
										, minBounds, maxBounds, edge_sol1[0], edge_sol2[0], quadraticFactors );
	if (m_debugOutput) 
	{
		cout << " #solutions: " << numSolutionsEdge[0] << endl; cout << endl;
	}

	edge.Set(dummy, maxBounds.GetY(), in_slicePositionZ);
	if (m_debugOutput) cout << "X_axis, top-edge: " << edge << " "; 
	numSolutionsEdge[1] = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, edge, Axis_X
										, minBounds, maxBounds, edge_sol1[1], edge_sol2[1], quadraticFactors );
	if (m_debugOutput) 
	{
		cout << " #solutions: " << numSolutionsEdge[1] << endl; cout << endl;
	}

	edge.Set(maxBounds.GetX(), dummy, in_slicePositionZ);
	if (m_debugOutput) cout << "Y_axis, right-edge: " << edge << " "; 
	numSolutionsEdge[2] = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, edge, Axis_Y
										, minBounds, maxBounds, edge_sol1[2], edge_sol2[2], quadraticFactors );
	if (m_debugOutput) 
	{
		cout << " #solutions: " << numSolutionsEdge[2] << endl; cout << endl;
	}

	edge.Set(dummy, minBounds.GetY(), in_slicePositionZ);
	if (m_debugOutput) cout << "X_axis, bottom-edge: " << edge << " ";
	numSolutionsEdge[3] = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, edge, Axis_X
										, minBounds, maxBounds, edge_sol1[3], edge_sol2[3], quadraticFactors );
	if (m_debugOutput) 
	{
		cout << " #solutions: " << numSolutionsEdge[3] << endl; cout << endl;
	}

	int totSolutions = numSolutionsEdge[0] + numSolutionsEdge[1] + numSolutionsEdge[2] + numSolutionsEdge[3];

	// /*
	if (m_debugOutput)
	{
		cout << endl;
		cout << "Looking for intersections with edges: [" << minBounds << " ; " << maxBounds << "] " << endl;  
		cout << "WILDERMAN march, total #solutions: " << totSolutions << endl;
		for (int ednr = 0; ednr < 4; ednr++)
		{
			cout << ", numSol" << ednr << ": " << numSolutionsEdge[ednr] 
				<< " with: " << edge_sol1[ednr] << " and: " << edge_sol2[ednr] << endl;
		}
		cout << endl;
	}
	// */

	// 1.b check that #number of solutions is even
	if (totSolutions % 2 != 0)
	{
		// TODO....
		// this, apparantly, can happen too (just on the border of the FOV)
		/*		
		cout << "ERROR! WildermanMarch(), #tot solutions: " << totSolutions << endl;
		for (int ednr = 0; ednr < 4; ednr++)
		{
			cout << ", numSol" << ednr << ": " << numSolutionsEdge[ednr] 
				 << " with: " << edge_sol1[ednr] << " and: " << edge_sol2[ednr] << endl;
		}
		exit(1);
		*/
		return scenario;
	}

	// 2. IdentifyIntersectionCase
	//		See copy below in file, maybe this is not necessary....
	// 2.1 Case 1 (special case)
	if (totSolutions == 0)
	{
		// check whether ellips is inside FOV
		// 2.1.a: intersection coneAxis with plane "in_slicePosition"
		double deltaZ = in_slicePositionZ - in_scatPos.GetZ();
		double factor = deltaZ / in_coneAxis.GetZ();
		C3Vector axisIntersection = in_scatPos + ( in_coneAxis * factor );
 		if ( !IsInsideBounds(minBounds, maxBounds, axisIntersection) )
		{
			// axis out of fov, so ellips too
			scenario = 0;
		}
		else
		{
			scenario = 1;

			if (doFillBins) 
			{
				// if the axis is inside, check if there are two intersection x points on y_slice, z_slice
				C3Vector sol1, sol2;
				int num = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, axisIntersection, Axis_X
												, minBounds, maxBounds, sol1, sol2, quadraticFactors );
				// assert( num == 2 );
				/*
				cout << "SCENARIO 1, find sols for axisIntersection: " << axisIntersection 
					<< " numSols: " << num << " SOL1: " << sol1 << " SOL2: " << sol2 << endl;
				*/
				/*
				// Just for cross-checking
				{
					C3Vector dum1, dum2;
					QuadraticFactors quadDum;
					C3Vector axisDum( 0.5*(sol1.GetX() + sol2.GetX()), 0, axisIntersection.GetZ());
					int dum = FindIntersectionsInFOV(in_scatPos, in_coneAxis, in_comptonAngle, axisDum, Axis_Y
												, minBounds, maxBounds, dum1, dum2, quadDum );
					cout << "SCENARIO 1, CROSS-CHECK, find sols for dummy axis Intersection: " << axisDum 
						<< " numSols: " << dum << " SOL1: " << dum1 << " SOL2: " << dum2 << endl;
				}
				*/

				if (num == 2 )
				{
					March( in_scatPos, in_coneAxis, in_comptonAngle, sol1, sol2
									, in_fieldOfView, io_fovidxVec, SEARCH_NEXT_X );
				}
			}
		}
	}
	else if (totSolutions == 2)
	{
		C3Vector sol1; 
		C3Vector sol2;
		bool foundpair = false;
		int firstedgenr = -1;
		for (int edgeidx = 0; edgeidx < 4 && !foundpair; edgeidx++)
		{
			// case 5: twice in one side, other sides zero
			if (numSolutionsEdge[edgeidx] == 2)
			{
				sol1 = edge_sol1[edgeidx];
				sol2 = edge_sol2[edgeidx];
				firstedgenr = edgeidx;
				foundpair = true;
				scenario = 5;
			}
			// case 2 or case "3": one in one side, and one in adjacent or opposite side
			else if (numSolutionsEdge[edgeidx] == 1)
			{
				sol1 = edge_sol1[edgeidx];				// the first solution from edge with index "edgeidx"
				firstedgenr = edgeidx;
				for (int edgeidx2 = edgeidx+1; edgeidx2 < 4 && !foundpair; edgeidx2++)
				{
					if (edgeidx2 != edgeidx && (numSolutionsEdge[edgeidx2] == 1))
					{
						sol2 = edge_sol1[edgeidx2];		// the first solution from edge with index "edgeidx2"
						foundpair = true;
						// 		scenario = (edgeidx2 = edgeidx+2) ? 3 : 2;
						// Wilderman calls the case of two opposite points "3" but 
						// I prefer to reserve "scenario 3" for the case of 4 total intersections
						scenario = 2;
					}
				}
			}	
		}
		assert( foundpair );
		assert( firstedgenr != -1);

		SEARCH_DIRECTION searchDirection;
		if ( firstedgenr == 0 )
			searchDirection = SEARCH_NEXT_X;
		else if ( firstedgenr == 1 )
			searchDirection = SEARCH_PREV_Y;
		else if ( firstedgenr == 2 )
			searchDirection = SEARCH_PREV_X;
		else if ( firstedgenr == 3 )
			searchDirection = SEARCH_NEXT_Y;
		else
		{
			cout << "ERROR! WildermanMarch(), #tot solutions: " << totSolutions << endl;
			for (int ednr = 0; ednr < 4; ednr++)
			{
				cout << ", numSol" << ednr << ": " << numSolutionsEdge[ednr] 
					 << " with: " << edge_sol1[ednr] << " and: " << edge_sol2[ednr] << endl;
				cout << "search Direction: " << searchDirection << endl;
			}
			assert( false );
		}

		// 3. March from first point to second point etc....
		if (doFillBins)
		{
			// cout << "SO MARCHING from sol1: " << sol1 << " to sol2: " << sol2 << endl;
			March( in_scatPos, in_coneAxis, in_comptonAngle, sol1, sol2, in_fieldOfView, io_fovidxVec, searchDirection );
		}
	}
	else if (totSolutions == 4)
	{
		// scenario 3 
		if (    ( numSolutionsEdge[0] == 2 && numSolutionsEdge[1] == 1 && numSolutionsEdge[3] == 1 )
			 || ( numSolutionsEdge[1] == 2 && numSolutionsEdge[0] == 1 && numSolutionsEdge[2] == 1 )
			 || ( numSolutionsEdge[2] == 2 && numSolutionsEdge[1] == 1 && numSolutionsEdge[3] == 1 )
			 || ( numSolutionsEdge[3] == 2 && numSolutionsEdge[0] == 1 && numSolutionsEdge[2] == 1 ) )
		{
			scenario = 3;
		}
		// scenario 4 
		else if (    ( numSolutionsEdge[0] + numSolutionsEdge[2] == 3 )
			      || ( numSolutionsEdge[1] + numSolutionsEdge[3] == 3 ) )
		{
			scenario = 4;
		}
		// scenario 6 
		else if (    ( numSolutionsEdge[0] == 2 && numSolutionsEdge[2] == 2 )
			      || ( numSolutionsEdge[1] == 2 && numSolutionsEdge[3] == 2 ) )
		{
			scenario = 6;
		}
		// scenario 10
		else if (    ( numSolutionsEdge[0] == 2 && numSolutionsEdge[1] == 2 )
			      || ( numSolutionsEdge[1] == 2 && numSolutionsEdge[2] == 2 )
			      || ( numSolutionsEdge[2] == 2 && numSolutionsEdge[3] == 2 )
			      || ( numSolutionsEdge[3] == 2 && numSolutionsEdge[0] == 2 ) )
		{
			scenario = 10;
		}
		else if (    ( numSolutionsEdge[0] == 1 && numSolutionsEdge[1] == 1 )
			      && ( numSolutionsEdge[2] == 1 && numSolutionsEdge[3] == 1 ) )
		{
			scenario = 11;
		}
		else
		{
			// /*
			cout << "ERROR! WildermanMarch(), #tot solutions: " << totSolutions << endl;
			for (int ednr = 0; ednr < 4; ednr++)
			{
				cout << ", numSol" << ednr << ": " << numSolutionsEdge[ednr] 
					 << " with: " << edge_sol1[ednr] << " and: " << edge_sol2[ednr] << endl;
			}
			//	assert( false );
			// */
		}
	}
	// totSolutions is 6 or 8
	else
	{
		if (totSolutions == 6)
			scenario = 7;
		else if (totSolutions == 8)
			scenario = 8;
		else
			scenario = 9;
	}

	return scenario;
}

// ****************************************************************************************************

void March( const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
			, const C3Vector& in_intersect1, const C3Vector& in_intersect2, const CVIPFieldOfView& in_fieldOfView
			, std::vector<unsigned int>& io_fovidxVec, SEARCH_DIRECTION in_searchDirection )
{
	// 1. Start up
    // begin point
	unsigned int fovidx = in_fieldOfView.GetVoxelIndex( in_intersect1 );
	io_fovidxVec.clear();
	io_fovidxVec.push_back( fovidx );

    // end-point
	unsigned int end_fovidx = in_fieldOfView.GetVoxelIndex( in_intersect2 );
	
	// /*
	if (m_debugOutput)
	{
		cout << "START pos: " << in_intersect1 << " with fovidx: " << fovidx 
			<< " END pos: " << in_intersect2 << " with fovidx: " << end_fovidx << endl;
	}
	// */
	if (fovidx == end_fovidx)
	{
		return;
	}

	// Call submarch
	SubMarch( fovidx, in_scatPos, in_coneAxis, in_coneAngle, in_intersect1, in_intersect2
					, in_fieldOfView, io_fovidxVec, in_searchDirection, end_fovidx );
}

// ****************************************************************************************************

void SubMarch( int in_fovidx, const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
			, const C3Vector& in_pnt1, const C3Vector& in_pnt2, const CVIPFieldOfView& in_fieldOfView
			, std::vector<unsigned int>& io_fovidxVec, SEARCH_DIRECTION in_searchDirection
			, unsigned int in_end_fovidx )
{

    int currentidxX, currentidxY, currentidxZ;
	in_fieldOfView.GetVoxelIndices( in_fovidx, currentidxX, currentidxY, currentidxZ );

    C3Vector edge1, edge2, sol1, sol2, sol3, sol4;
    int uppBin, lowBin;
	AXIS axisSearch = (in_searchDirection == SEARCH_NEXT_X) || (in_searchDirection == SEARCH_PREV_X) ? Axis_X : Axis_Y;
	double unknown = 0.0, currentPosCoeff;

	const C3Vector& minBounds = in_fieldOfView.GetLowerEdge();
	const C3Vector& maxBounds = in_fieldOfView.GetUpperEdge();
	
    if ( axisSearch == Axis_X )	   // SEARCH_NEXT_X || SEARCH_PREV_X  = search for X solutions on upper and lower Y border
    {
		// A search on the X-axis means: 
		// we look for new values of X, along two (upper/lower) horizontal lines with known Y coefficients ("upp" and "low").
		uppBin = currentidxY + 1;
		double yUpp = minBounds.GetY() + uppBin * in_fieldOfView.GetVoxelSize().GetY();
		edge1.Set( unknown, yUpp, in_pnt1.GetZ() );

		lowBin = currentidxY;
		double yLow = minBounds.GetY() + lowBin * in_fieldOfView.GetVoxelSize().GetY();
		edge2.Set( unknown, yLow, in_pnt1.GetZ() );

		currentPosCoeff = in_pnt1.GetX();
    }
	else							// SEARCH_NEXT_Y || SEARCH_PREV_Y	-> search for Y solutions on left and right X border
	{
		// A search on the Y-axis means: 
		// we look for new values of Y, along two (left/right) vertical lines with known X coefficients ("upp" and "low").
		uppBin = currentidxX + 1;
		double xUpp = minBounds.GetX() + uppBin * in_fieldOfView.GetVoxelSize().GetX();
		edge1.Set( xUpp, unknown, in_pnt1.GetZ() );

		lowBin = currentidxX;
		double xLow = minBounds.GetX() + lowBin * in_fieldOfView.GetVoxelSize().GetX();
		edge2.Set( xLow, unknown, in_pnt1.GetZ() );

		currentPosCoeff = in_pnt1.GetY();
	}

	/*
	cout << " ------------------------------------------------- " << endl;
	cout << "subMarch, search Dir: ";
	PrintSearchDirection( in_searchDirection );
	cout << " currentPosCoeff: " << currentPosCoeff << endl;
	*/

	// Find intersections of the Compton cone with the edges (i.e. the horizontal/vertical lines as defined before)
	QuadraticFactors quadraticFactors1, quadraticFactors2;
	int numSolutions1 = FindIntersectionsInFOV( in_scatPos, in_coneAxis, in_coneAngle, edge1, axisSearch
						, minBounds, maxBounds, sol1, sol2, quadraticFactors1 );
	int numSolutions2 = FindIntersectionsInFOV( in_scatPos, in_coneAxis, in_coneAngle, edge2, axisSearch
						, minBounds, maxBounds, sol3, sol4, quadraticFactors2 );

	// Skip solutions that are same as current position
	CheckEquals(numSolutions1, in_pnt1, sol1, sol2);
	CheckEquals(numSolutions2, in_pnt1, sol3, sol4);

	/*
	cout << "(after check equals) for edge1: " << edge1 << " numSols: " 
		 << numSolutions1 << " SOL1: " << sol1 << " SOL2: " << sol2 << endl;
	cout << "      		      and for edge2: " << edge2 << " numSols: " 
		 << numSolutions2 << " SOL3: " << sol3 << " SOL4: " << sol4 << endl;	
	*/

	// Find which solution is the closest point to the current point
    int whichSolution = 0;
	if (numSolutions1 + numSolutions2 > 0)
	{
		whichSolution = FindNearestPoint(in_searchDirection, currentPosCoeff
								, numSolutions1, sol1, sol2, numSolutions2, sol3, sol4);
	}
	
	//	/*
	if (m_debugOutput)
	{
		cout << "Number of solutions: " << numSolutions1 << " + " << numSolutions2 << endl;
		cout << "Nearest point sol: " << whichSolution << "   ";
		cout << "SO, the solution is: ";
		if (whichSolution == 0)
			cout << "FINAL point: " << in_pnt2 << endl;
		else if (whichSolution == 1)
			cout << sol1 << endl;
		else if (whichSolution == 2)
			cout << sol2 << endl;
		else if (whichSolution == 3)
			cout << sol3 << endl;
		else if (whichSolution == 4)
			cout << sol4 << endl;
	}
	// 	*/

	// Get solution coefficients and prepare variabeles for the next sub march
	//
	C3Vector nextPnt;	
			// This is the final solution, corresponding to the next intersection in the FOV starting from the current position
			// All the pixels between the current position and the "next position" will be fillled
	int nextidxX, nextidxY, dummyidx, nextidxQ;
	SEARCH_DIRECTION searchDirection = in_searchDirection;
	if (whichSolution == 0)
	{
		// In this case, we only going to fill all pixels from pnt1 to pnt2 and not calling another march.
		in_fieldOfView.GetVoxelIndices(in_pnt2, nextidxX, nextidxY, dummyidx);
		nextidxQ = ( axisSearch == Axis_X ) ? nextidxX : nextidxY;
		searchDirection = SEARCH_NONE;
	}
	// Found a solution on a border with larger X-or-Y coefficient than the current point
	else if (whichSolution == 1 || whichSolution == 2)
	{
		if (whichSolution == 1)
		{	
			nextPnt = sol1;
			in_fieldOfView.GetVoxelIndices(sol1, nextidxX, nextidxY, dummyidx);
		}
		else if (whichSolution == 2)
		{
			nextPnt = sol2;
			in_fieldOfView.GetVoxelIndices(sol2, nextidxX, nextidxY, dummyidx);
		}

		if ( axisSearch == Axis_X )
		{
			if ( nextidxY != uppBin && uppBin != in_fieldOfView.GetNumVoxelsY() )
			{
				cout << "SOMETHING WRONG, nextidxY: " << nextidxY << " uppBin: " << uppBin
					 << " fov.numVoxels: " << in_fieldOfView.GetNumVoxelsY() << endl;
				cout << "MOST LIKELY THIS IS A ROUNDING OFF MISTAKE (FOV voxels)" << endl;
				return;
			}
			// assert( nextidxY == uppBin || uppBin == in_fieldOfView.GetNumVoxelsY() );
			nextidxQ = nextidxX;
		}
		else
		{
			if ( nextidxX != uppBin && uppBin != in_fieldOfView.GetNumVoxelsX() )
			{
				cout << "SOMETHING WRONG, nextidxX: " << nextidxX << " uppBin: " << uppBin
					 << " fov.numVoxels: " << in_fieldOfView.GetNumVoxelsX() << endl;
				cout << "MOST LIKELY THIS IS A ROUNDING OFF MISTAKE (FOV voxels)" << endl;
				return;
			}
			// assert( nextidxX == uppBin || uppBin == in_fieldOfView.GetNumVoxelsX() );
			nextidxQ = nextidxY;
		}

		// If we come from lower X or Y, than the next search will be for an even higher value
		// 		unless (rare case) we just touch the edge but do not cross it...
		double discriminant = (quadraticFactors1.GetbFactor() * quadraticFactors1.GetbFactor()) 
					 	- (4 * quadraticFactors1.GetaFactor() * quadraticFactors1.GetcFactor());
		if (discriminant != 0) // Otherwise, DIS == 0: only 1 solution, i.e., not crossing but only touching this line
		{
			searchDirection = (axisSearch == Axis_X) ? SEARCH_NEXT_Y : SEARCH_NEXT_X;
		}
	}
	// Found a solution on a border with smaller X-or-Y coefficient than the current point
	else if (whichSolution == 3 || whichSolution == 4)
	{
		if (whichSolution == 3)
		{
			nextPnt = sol3;
			in_fieldOfView.GetVoxelIndices(sol3, nextidxX, nextidxY, dummyidx);
		}
		else if (whichSolution == 4)
		{
			nextPnt = sol4;
			in_fieldOfView.GetVoxelIndices(sol4, nextidxX, nextidxY, dummyidx);
		}

		if ( axisSearch == Axis_X )
		{

			if (nextidxY != lowBin)
			{
				cout << "AxisX, nextidxY: " << nextidxY << " lowBin: " << lowBin << endl;
				cout << "MOST LIKELY THIS IS A ROUNDING OFF MISTAKE (FOV voxels)" << endl;
				return;
				// nextidxY = lowBin;
			}

			assert(nextidxY == lowBin);
			nextidxQ = nextidxX;
			nextidxY--;
			if (nextidxY < 0) nextidxY = 0;
		}
		else
		{

			if (nextidxX != lowBin)
			{
				cout << "AxisY, nextidxX: " << nextidxX << " lowBin: " << lowBin << endl;
				cout << "MOST LIKELY THIS IS A ROUNDING OFF MISTAKE (FOV voxels)" << endl;
				return;
				// nextidxX = lowBin;
			}

			assert(nextidxX == lowBin);
			nextidxQ = nextidxY;
			nextidxX--;
			if (nextidxX < 0) nextidxX = 0;
		}

		// If we come from higher X or Y, than the next search will be for an even lower value
		// unless (rare case) we just touch the edge but do not cross it...
		double discriminant = (quadraticFactors2.GetbFactor() * quadraticFactors2.GetbFactor()) 
						- (4 * quadraticFactors2.GetaFactor() * quadraticFactors2.GetcFactor());
		if (discriminant != 0)	// Otherwise, DIS == 0: only 1 solution, i.e., not crossing but only touching this line
		{
			searchDirection = (axisSearch == Axis_X) ? SEARCH_PREV_Y : SEARCH_PREV_X;
		}
	}

    // Fill bins
	bool doStop = false;
	FillBins( axisSearch, currentidxX, currentidxY, currentidxZ, nextidxQ
							, in_fieldOfView, io_fovidxVec, doStop, in_end_fovidx );

	if (m_debugOutput)
	{
		cout << "searchDirection: " << searchDirection << " =? NONE: " << (searchDirection == SEARCH_NONE)
		     << " doStop?: " << doStop << endl;
	}
	
    // Call next sub-march
	if (searchDirection != SEARCH_NONE && !doStop )
	{
		int fovidx;		// the FOV index of the next starting position
		// When we searched for a new position on the X-axis, the Y and Z coefficients are equal to previous position
		if (axisSearch == Axis_X)
		{
			fovidx = in_fieldOfView.GetVoxelIndex( nextidxQ, currentidxY, currentidxZ );
		}
		// When we searched for a new position on the Y-axis, the X and Z coefficients are equal to previous position
		else if (axisSearch == Axis_Y) 
		{
			fovidx = in_fieldOfView.GetVoxelIndex( currentidxX, nextidxQ, currentidxZ );
		}

		SubMarch( fovidx, in_scatPos, in_coneAxis, in_coneAngle, nextPnt, in_pnt2
						, in_fieldOfView, io_fovidxVec, searchDirection, in_end_fovidx );
	}
}

// ****************************************************************************************************

// Loop from current current bin to nearest bin

void
FillBins( AXIS in_axis, int in_currentidxX, int in_currentidxY, int in_currentidxZ, int in_nextidxQ
         , const CVIPFieldOfView& in_fieldOfView, std::vector<unsigned int>& io_fovidxVec
         , bool& io_doStop, unsigned int in_end_fovidx )
{
	io_doStop = false;
    int currentidxQ = (in_axis == Axis_X) ? in_currentidxX : in_currentidxY;
	if (currentidxQ != in_nextidxQ)
	{
		unsigned int fovidx;
		std::vector<unsigned int>::const_iterator iter;
		unsigned int startidx = (in_nextidxQ > currentidxQ) ? currentidxQ+1 : in_nextidxQ;
		unsigned int endidx = (in_nextidxQ > currentidxQ) ? in_nextidxQ : currentidxQ-1;
	
		for (int iq = startidx; iq <= endidx; iq++)
		{
			if (in_axis == Axis_X)
			{
				fovidx = in_fieldOfView.GetVoxelIndex( iq, in_currentidxY, in_currentidxZ );
			}
			else if (in_axis == Axis_Y)
			{
				fovidx = in_fieldOfView.GetVoxelIndex( in_currentidxX, iq, in_currentidxZ );
			}
	
			iter = std::find(io_fovidxVec.begin(), io_fovidxVec.end(), fovidx);
			if (iter == io_fovidxVec.end()) // if NOT found already....
			{
				io_fovidxVec.push_back( fovidx );
			}
			else
			{
				// If the final position has been found, STOP
				// if ( std::find(io_fovidxVec.begin(), io_fovidxVec.end(), in_end_fovidx) != io_fovidxVec.end() )
				//	
				// We cannot apply this condition. Even when the final position has NOT been found, we have to stop
				// 	because we are stuck and will end up in an endless loop...
				{
					io_doStop = true;
					if (m_debugOutput)
					{
						cout << " STOP, because solution already found before; fovidx: " << fovidx << endl;
					}
				}
			}
		}
	}
	// If not the first search (which sometimes ends up in same bin...)
	// TODO: Does this happen? Why?
	else if (io_fovidxVec.size() > 1)	
	{
		if (m_debugOutput)
		{		
			cout << " STOP, because io_fovidxVec.size(): " << io_fovidxVec.size() 
				 << " and currentidxQ == in_nextidxQ: " << currentidxQ << "  " << in_nextidxQ << endl;
		}
		io_doStop = true;
	}
}

// ****************************************************************************************************

int
FindNearestPoint( SEARCH_DIRECTION in_searchDirection, const double& in_currentPos, 
		int in_numSolutions1, const C3Vector& in_sol1, const C3Vector& in_sol2, 
		int in_numSolutions2, const C3Vector& in_sol3, const C3Vector& in_sol4)
{
	// Simpler variables
	double sol1Pos = (in_searchDirection == SEARCH_NEXT_X || in_searchDirection == SEARCH_PREV_X) ? in_sol1.GetX() : in_sol1.GetY();
	double sol2Pos = (in_searchDirection == SEARCH_NEXT_X || in_searchDirection == SEARCH_PREV_X) ? in_sol2.GetX() : in_sol2.GetY();
	double sol3Pos = (in_searchDirection == SEARCH_NEXT_X || in_searchDirection == SEARCH_PREV_X) ? in_sol3.GetX() : in_sol3.GetY();
	double sol4Pos = (in_searchDirection == SEARCH_NEXT_X || in_searchDirection == SEARCH_PREV_X) ? in_sol4.GetX() : in_sol4.GetY();

	// cout << "sols, in FindNearestPos, 1: " << sol1Pos << " " << sol2Pos << " " << sol3Pos << " " << sol4Pos << endl;

	Comparator::COMPARISON comparisonCurrentToNext = (in_searchDirection == SEARCH_NEXT_X || in_searchDirection == SEARCH_NEXT_Y) 
				? Comparator::COMPARE_SMALLER : Comparator::COMPARE_BIGGER;
	Comparator comparator( comparisonCurrentToNext ); 

	int sol12 = 0; 
	int sol34 = 0;
	double sol12Pos, sol34Pos, solPos;

	if (in_numSolutions1 == 2)
	{
		sol12 = CompareThreePoints( in_currentPos, sol1Pos, sol2Pos, comparator );
	}
	else if (in_numSolutions1 == 1)
	{
		sol12 = comparator.Compare(in_currentPos, sol1Pos) ? 1 : 0;
	}

	if (sol12 != 0)
		sol12Pos = (sol12 == 1) ? sol1Pos : sol2Pos;

	if (in_numSolutions2 == 2)
	{
		sol34 = CompareThreePoints( in_currentPos, sol3Pos, sol4Pos, comparator );
	}
	else if (in_numSolutions2 == 1)
	{
		sol34 = comparator.Compare(in_currentPos, sol3Pos) ? 1 : 0;
	}

	if (sol34 != 0)
	{
		sol34 += 2;	// values are 3 or 4 (instead of 1 or 2)
		sol34Pos = (sol34 == 3) ? sol3Pos : sol4Pos;
	}

	// if no more crossings found, we might be at the last stretch (i.e. from here to end-point)
	if (sol12 == 0 && sol34 == 0)
	{
		return 0;
	}

	int sol = 0;
	if (sol12 != 0 && sol34 != 0)
	{
		sol = CompareThreePoints( in_currentPos, sol12Pos, sol34Pos, comparator );
		if (sol == 1) 
			sol = sol12;
		else if (sol == 2) 
			sol = sol34;
	}
	else 
	{
		sol = (sol12 != 0) ? sol12 : sol34;
	}
	return sol;
}

// ****************************************************************************************************

void
CheckEquals(int& io_numSolutions, const C3Vector& in_currentPos, C3Vector& io_sol1, C3Vector& io_sol2)
{
	if (io_numSolutions == 2)
	{
		if (io_sol2 == in_currentPos)
		{
			io_numSolutions = 1;
			if (io_sol1 == in_currentPos)
			{
				io_numSolutions = 0;
			}
		}
		else if (io_sol1 == in_currentPos)
		{
			io_numSolutions = 1;
			io_sol1 = io_sol2;
		}
	}
	else if (io_numSolutions == 1 && io_sol1 == in_currentPos)
	{
		io_numSolutions = 0;
	}
}

// ****************************************************************************************************

int
CompareThreePoints( const double& in_currentPos, const double& in_sol1Pos, const double& in_sol2Pos, const Comparator& in_comparator )
{
	if ( !in_comparator.Compare(in_currentPos, in_sol1Pos) && !in_comparator.Compare(in_currentPos, in_sol2Pos) )
	{
		return 0;
	}
	else if ( in_comparator.Compare(in_currentPos, in_sol1Pos) && !in_comparator.Compare(in_currentPos, in_sol2Pos) )
	{
		return 1;
	}
	else if ( !in_comparator.Compare(in_currentPos, in_sol1Pos) && in_comparator.Compare(in_currentPos, in_sol2Pos) )
	{
		return 2;
	}
	else if ( in_comparator.Compare(in_currentPos, in_sol1Pos) && in_comparator.Compare(in_currentPos, in_sol2Pos) )
	{	
		return (in_comparator.Compare(in_sol1Pos, in_sol2Pos)) ? 1 : 2;
	}
	else
		assert( false );
}

// ****************************************************************************************************

/*
	Scott J. Wilderman: 
		equation (1) from "Fast Algorithm for List Mode Back-Projection of Compton Scatter Camera Data"
	IEEE Transactions on Nucl. Science Vol. 45 1998
*/
int FindIntersectionsInFOV(const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_coneAngle
        , const C3Vector& in_edge, AXIS in_toSolve, const C3Vector& in_minBounds, const C3Vector& in_maxBounds
		, C3Vector& io_intersect1, C3Vector& io_intersect2, QuadraticFactors& io_quadraticFactors )
{
	// Call function to calculate the a,b,c factors for the quadratic equation....
	FromConeEquationToQuadraticFactor( in_coneAxis, in_scatPos, in_coneAngle, in_toSolve, in_edge, io_quadraticFactors );
	if (m_debugOutput) 
	{
		cout << "FindIntersectionsInFOV, factors: " << io_quadraticFactors.GetaFactor() << " " << io_quadraticFactors.GetbFactor() 
             << " " << io_quadraticFactors.GetcFactor() << endl;
	}

    double q1 = 0.0, q2 = 0.0;
    int numSolutionsQ = solveQuadraticEquation(io_quadraticFactors, q1, q2);
	if (m_debugOutput) 
	{
		cout << "FindIntersectionsInFOV, numSolutionsQ: " << numSolutionsQ << " q1: " << q1 << " q2: " << q2 << endl;
	}

    int numSolutions = 0;
    io_intersect1.Set(-1,-1,-1);
    io_intersect2.Set(-1,-1,-1);

    if (numSolutionsQ > 0)
    {
		// Wildermann: 
		// solution to quadratic equation q1 = (x - x1, y - y1, z - z1); with (x1, y1, z1) = 1st hit
		// So, position of the solution is (x, y, z) = q1 + scatPos
		double x = (in_toSolve == Axis_X) ? q1 + in_scatPos.GetX() : in_edge.GetX();
        double y = (in_toSolve == Axis_Y) ? q1 + in_scatPos.GetY() : in_edge.GetY();
        double z = (in_toSolve == Axis_Z) ? q1 + in_scatPos.GetZ() : in_edge.GetZ();
		C3Vector tmpSolPos(x, y, z);

		C3Vector direction = tmpSolPos - in_scatPos;
		double inprod = direction.GetScalarProductCosAngle( in_coneAxis );
		bool cosPositive = (cos(in_coneAngle) > 0);
		bool inprodPositive = (inprod > 0);
		bool isOK = (cosPositive && inprodPositive) || (!cosPositive && !inprodPositive);

		if (m_debugOutput)
		{
			cout << "Found solution, with position: " << tmpSolPos << " isOK: " << isOK 
				<< " inFOV: " << IsInsideBounds(in_minBounds, in_maxBounds, tmpSolPos) << endl;
		}
		if ( IsInsideBounds(in_minBounds, in_maxBounds, tmpSolPos) && isOK )
		{
			io_intersect1 = tmpSolPos;
        	numSolutions++;
		}
        if (numSolutionsQ > 1)
        {
			x = (in_toSolve == Axis_X) ? q2 + in_scatPos.GetX() : in_edge.GetX();
			y = (in_toSolve == Axis_Y) ? q2 + in_scatPos.GetY() : in_edge.GetY();
			z = (in_toSolve == Axis_Z) ? q2 + in_scatPos.GetZ() : in_edge.GetZ();
			tmpSolPos.Set(x, y, z);

			direction = tmpSolPos - in_scatPos;
			inprod = direction.GetScalarProductCosAngle( in_coneAxis );
			inprodPositive = (inprod > 0);
			isOK = (cosPositive && inprodPositive) || (!cosPositive && !inprodPositive);

			if (m_debugOutput)
			{
				cout << "Found other solution, with position: " << tmpSolPos << " isOK: " << isOK 
					<< " inFOV: " << IsInsideBounds(in_minBounds, in_maxBounds, tmpSolPos) << endl;
			}			
			
			if ( IsInsideBounds(in_minBounds, in_maxBounds, tmpSolPos) && isOK )
			{
				if (numSolutions == 0)
					io_intersect1 = tmpSolPos;
				else
        	    	io_intersect2 = tmpSolPos;
            	numSolutions++;
			}
        }

    }
	return numSolutions;
}

// ****************************************************************************************************

/*
	Scott J. Wilderman: 
		equation (1) from "Fast Algorithm for List Mode Back-Projection of Compton Scatter Camera Data"
	IEEE Transactions on Nucl. Science Vol. 45 1998
*/
void 
FromConeEquationToQuadraticFactor( const C3Vector& in_coneAxis, const C3Vector& in_scatPos, const double& in_coneAngle 
	, AXIS in_toSolve, const C3Vector& in_edge, QuadraticFactors& io_quadraticFactors )
{
	double lambda = cos(in_coneAngle);
	double lambda2 = lambda*lambda;
	double inverselambda2 = 1.0/lambda2;

	// cout << "lambda: " << lambda << ", lambda^2: " << lambda2 << endl;

	// Make sure cone axis is normalized to length 1!!!!
	double len = in_coneAxis.GetLength();
	double factor = 1.0;
	if (len != 1.0)
	{
		factor = len > 0 ? (1.0 / len) : 0.0;
	}
    double nx = factor * in_coneAxis.GetX();
    double ny = factor * in_coneAxis.GetY();
    double nz = factor * in_coneAxis.GetZ();
	
	if (m_debugOutput)
	{	
		cout << "Cone axis: " << nx << " " << ny << " " << nz << " lambda2: " << lambda2 << endl;
	}
	
	double x1 = in_scatPos.GetX();
	double y1 = in_scatPos.GetY();
	double z1 = in_scatPos.GetZ();

	double T, R;
	if (in_toSolve == Axis_X) // solve X, so Y (and Z) is known
	{
		double y = in_edge.GetY();
		double zs = in_edge.GetZ();
		T = ny * (y - y1) + nz * (zs - z1);
		R = (y - y1)*(y - y1) + (zs - z1)*(zs - z1);
		//	testFnn = nx*nx; testFl2 = lambda2; testFnnMINFl2 = nx*nx - lambda2;
		//	io_quadraticFactors.SetaFactor( nx*nx - lambda2 );
		io_quadraticFactors.SetaFactor( doubleEquals(nx*nx, lambda2, 0.0001) ? 0.0 : nx*nx - lambda2 );
		io_quadraticFactors.SetbFactor( 2 * T * nx );
	}
	else if (in_toSolve == Axis_Y) // solve Y, so X (and Z) is known
	{
		double x = in_edge.GetX();
		double zs = in_edge.GetZ();
		T = nx * (x - x1) + nz * (zs - z1);
		R = (x - x1)*(x - x1) + (zs - z1)*(zs - z1);
		//	io_quadraticFactors.SetaFactor( ny*ny - lambda2 );
		io_quadraticFactors.SetaFactor( doubleEquals(ny*ny, lambda2, 0.0001) ? 0.0 : ny*ny - lambda2 );
		io_quadraticFactors.SetbFactor( 2 * T * ny );
	}
	else if (in_toSolve == Axis_Z) // solve Z, so X and Y are known
	{
		double x = in_edge.GetX();
		double y = in_edge.GetY();
		T = nx * (x - x1) + ny * (y - y1);
		R = (x - x1)*(x - x1) + (y - y1)*(y - y1);
		//	io_quadraticFactors.SetaFactor( nz*nz - lambda2 );
		io_quadraticFactors.SetaFactor( doubleEquals(nz*nz, lambda2, 0.0001) ? 0.0 : nz*nz - lambda2 );
		io_quadraticFactors.SetbFactor( 2 * T * nz );
	}

	double T2 = T*T;
	double lambda2R = lambda2 * R;

	io_quadraticFactors.SetcFactor( doubleEquals(T2, lambda2R, 0.0001) ? 0.0 : T2 - lambda2R );		// T*T - lambda2 * R;

	if (m_debugOutput)
	{
		cout << "T: " << T << ", T*T: " << T2 << ", R: " << R 
             << ", lambda2 * R: " << lambda2R << " (T^2 - l^2 R:) " 
             << T2 - lambda2R << " diff: " << endl;
		cout << "Quadratic factors: " 
				 << io_quadraticFactors.GetaFactor() << "  "
				 << io_quadraticFactors.GetbFactor() << "  "
				 << io_quadraticFactors.GetcFactor() << "  " << endl;
		cout << "	facA, nx*nx: " << nx*nx << " ny*ny: " << ny*ny << " nz*nz: " << nz*nz << endl;
		cout << "	facB, 2*T*nx: " << 2*T*nx << " 2*T*ny: " << 2*T*ny << " 2*T*nz: " << 2*T*nz << endl;
		cout << "	facC, T2: " << T2 << " lambda2R: " << lambda2R << endl;
	}
}

// ****************************************************************************************************
/*
void
FromConeEquationToQuadraticFactorOld( const C3Vector& in_coneAxis, const C3Vector& in_scatPos, const double& in_coneAngle 
	, AXIS in_toSolve, const C3Vector& in_edge, QuadraticFactors& io_quadraticFactors )
{
	// Make sure cone axis is normalized to length 1!!!!
	double len = in_coneAxis.GetLength();
	double factor = 1.0;
	if (len != 1.0)
	{
		factor = len > 0 ? (1.0 / len) : 0.0;
	}
    double nx = factor * in_coneAxis.GetX();
    double ny = factor * in_coneAxis.GetY();
    double nz = factor * in_coneAxis.GetZ();

    double LTERM = 0.0;   // left hand known term
    double RTERM = 0.0;   // right hand known term

    if (in_toSolve == Axis_Y || in_toSolve == Axis_Z)    	// calculate unknown y or z from known x
    {
        LTERM = nx * (in_edge.GetX() - in_scatPos.GetX());
        RTERM = (in_edge.GetX() - in_scatPos.GetX()) * (in_edge.GetX() - in_scatPos.GetX());
    }
    if (in_toSolve == Axis_X || in_toSolve == Axis_Z)       // calculate unknown x or z from known y
    {
        LTERM += ny * (in_edge.GetY() - in_scatPos.GetY());
        RTERM += (in_edge.GetY() - in_scatPos.GetY()) * (in_edge.GetY() - in_scatPos.GetY());
    }
    if (in_toSolve == Axis_X || in_toSolve == Axis_Y)       // calculate unknown x or y from known z
	{
    	LTERM += nz * (in_edge.GetZ() - in_scatPos.GetZ());
    	RTERM += (in_edge.GetZ() - in_scatPos.GetZ()) * (in_edge.GetZ() - in_scatPos.GetZ());
	}

	// cout << "LTERM: " << LTERM << ", RTERM: " << RTERM << endl;

    // root-square formula
	double cosAngle = cos(in_coneAngle);
    double cosSQ = cosAngle * cosAngle;
	io_quadraticFactors.SetcFactor( LTERM*LTERM - cosSQ * RTERM );
    if (in_toSolve == Axis_Y)
    {
        io_quadraticFactors.SetbFactor( 2 * LTERM * ny );
        io_quadraticFactors.SetaFactor( ny*ny - cosSQ );
    }
    else if (in_toSolve == Axis_X)
    {
        io_quadraticFactors.SetbFactor( 2 * LTERM * nx );
        io_quadraticFactors.SetaFactor( nx*nx - cosSQ );
    }
    else if (in_toSolve == Axis_Z)
    {
        io_quadraticFactors.SetbFactor( 2 * LTERM * nz );
        io_quadraticFactors.SetaFactor( nz*nz - cosSQ );
    }
}
*/
// ****************************************************************************************************

void
ShowFovVector( const CVIPFieldOfView& in_fieldOfView, const std::vector<unsigned int>& in_fovidxVec )
{
	cout << "SHOW, fov idx list size: " << in_fovidxVec.size() << " with fovidxs: " << endl;
	for (int i = 0; i < in_fovidxVec.size(); i++)
	{
		cout << in_fovidxVec[i] << "   ";
	}
	cout << endl;

	cout << "DRAWING: " << endl;
	if ( in_fovidxVec.size() > 0 )
	{
		std::vector<unsigned int>::const_iterator iter;
		for (int iz = 0; iz < in_fieldOfView.GetNumVoxelsZ(); iz++)
		{
			for (int iy = in_fieldOfView.GetNumVoxelsY() - 1; iy >= 0; iy--)
			{
				for (int ix = 0; ix < in_fieldOfView.GetNumVoxelsX(); ix++)
				{
					cout << "------";
				}
				cout << endl;
				for (int ix = 0; ix < in_fieldOfView.GetNumVoxelsX(); ix++)
				{
					unsigned int currIdx = in_fieldOfView.GetVoxelIndex( ix, iy, iz );
					// cout << "[" << ix << ", " << iy << "] " << "currIdx: " << currIdx << endl;
		
					iter = std::find(in_fovidxVec.begin(), in_fovidxVec.end(), currIdx);
					if (iter != in_fovidxVec.end()) 	// if found....
						cout << "|  x  ";
					else								// if NOT found....
						cout << "|  .  ";
				}
				cout << "|" << endl;
			}
			for (int ix = 0; ix < in_fieldOfView.GetNumVoxelsX(); ix++)
			{
				cout << "------";
			}
			cout << endl;
		}
	}

	cout << endl;
}

// ****************************************************************************************************

void
CalculateEllipsFactors( const double& in_slicePositionZ   
			, const C3Vector& in_scatPos, const C3Vector& in_coneAxis, const double& in_comptonAngle
			, const CVIPFieldOfView& in_fieldOfView, EllipsParameters& io_ellipsParameters )
{
	double lambda = cos(in_comptonAngle);
	double lambda2 = lambda*lambda;
	double inverselambda2 = 1.0/lambda2;

	double len = in_coneAxis.GetLength();
	double factor = 1.0;
	if (len != 1.0)
	{
		factor = len > 0 ? (1.0 / len) : 0.0;
	}
    double nx = factor * in_coneAxis.GetX();
    double ny = factor * in_coneAxis.GetY();
    double nz = factor * in_coneAxis.GetZ();

	double x1 = in_scatPos.GetX();
	double y1 = in_scatPos.GetY();
	double z1 = in_scatPos.GetZ();

	double nx2minlam2 = nx*nx - lambda2;
	double ny2minlam2 = ny*ny - lambda2;
	double nz2minlam2 = nz*nz - lambda2;
	double z = (in_slicePositionZ - z1);

	io_ellipsParameters.A = nx2minlam2;
	io_ellipsParameters.B = ny2minlam2;
	io_ellipsParameters.C = nx*ny;
	io_ellipsParameters.D = 2 * (-1*nx2minlam2*x1 - nx*ny*y1 + nx*nz*z);
	io_ellipsParameters.E = 2 * (-1*ny2minlam2*y1 - nx*ny*x1 + ny*nz*z);
	io_ellipsParameters.F = nx2minlam2*x1*x1 + ny2minlam2*y1*y1 + 2*nx*ny*x1*y1 - 2*nx*nz*x1*z - 2*ny*nz*y1*z + nz2minlam2*z*z;

	// yeah, but maybe not...
	io_ellipsParameters.A = -1 * io_ellipsParameters.A * inverselambda2;
	io_ellipsParameters.B = -1 * io_ellipsParameters.B * inverselambda2;
	io_ellipsParameters.C = -1 * io_ellipsParameters.C * inverselambda2;
	io_ellipsParameters.D = -1 * io_ellipsParameters.D * inverselambda2;
	io_ellipsParameters.E = -1 * io_ellipsParameters.E * inverselambda2;
	io_ellipsParameters.F = -1 * io_ellipsParameters.F * inverselambda2;

	AdditionalParameters( io_ellipsParameters );
}

// ****************************************************************************************************

/*
Source: http://www.purplemath.com/modules/sqrellps.htm
*/

void
AdditionalParameters( EllipsParameters& io_ellipsParameters )
{
    io_ellipsParameters.axis_a = 0.0;
    io_ellipsParameters.axis_b = 0.0;

    io_ellipsParameters.centre_x = 0.0;
    io_ellipsParameters.centre_y = 0.0;

    if ( io_ellipsParameters.A != 0 )
    {
        io_ellipsParameters.centre_x = -0.5 * io_ellipsParameters.D / io_ellipsParameters.A;
    }
    if ( io_ellipsParameters.B != 0 )
    {
        io_ellipsParameters.centre_y = -0.5 * io_ellipsParameters.E / io_ellipsParameters.B;
    }

    if ( io_ellipsParameters.A != 0 && io_ellipsParameters.B != 0 )
    {
        double G = - io_ellipsParameters.F / (io_ellipsParameters.A * io_ellipsParameters.B);
        G += ( io_ellipsParameters.centre_x * io_ellipsParameters.centre_x ) / io_ellipsParameters.B;
        G += ( io_ellipsParameters.centre_y * io_ellipsParameters.centre_y ) / io_ellipsParameters.A;

        if (io_ellipsParameters.B * G > 0)
        {
            io_ellipsParameters.axis_a = sqrt(io_ellipsParameters.B * G);
        }
        if (io_ellipsParameters.A * G > 0)
        {
            io_ellipsParameters.axis_b = sqrt(io_ellipsParameters.A * G);
        }
    }
}












