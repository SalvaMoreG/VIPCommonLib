

#include "CVIPMapData.h"
#include "CVIPUtils.h"
#include "VIPconstants.h"
#include "CVIPImage.h"
#include "CVIPFieldOfView.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <algorithm>

enum PLANE_AXIS
{
	PLANE_AXIS_NONE = 0,
	PLANE_AXIS_X,
	PLANE_AXIS_Y,
	PLANE_AXIS_Z
};

using namespace std;

void
GetLorPlaneAxis( int in_firstidx, int in_lastidx, const CVIPFieldOfView& in_fieldOfView, PLANE_AXIS& out_planeAxis );

void
MakeFakeFovMap( const CVIPFieldOfView& in_fieldOfView, CVIPMapData<unsigned int>& io_eventFovDataMap );

// Parameters
int m_eventMapSize = 0;		// TODO: input
int m_numMapFiles = 0;		// TODO: input
double m_tubeRadius = 0.0; 	// TODO: input

bool m_makeFakeFovMap( false );
bool m_isIntersection( false );

void 
ReadParameters();

// ====================================================================================

int main( int argc, char *argv[] )
{
	CVIPFieldOfView fieldOfView;
	fieldOfView.Initialize( "tor_fov_parameters.conf" );

	if (argc > 1)
	{
		string str = argv[1];
cout << "input: " << str << endl;
		if (str == "-makefakefov")
		{
			m_makeFakeFovMap = true;
		}
		else if (str == "-intersection")
		{
			m_isIntersection = true;
		}
		else
		{
			cout << "MakeTubeOfResponseFovMaps -makefakefov" << endl;
			cout << "MakeTubeOfResponseFovMaps -intersection" << endl;
			return 0;
		}
	}

	// read TOR parameters
	ReadParameters();

	// FOV event Maps...
	CVIPMapData<unsigned int> eventFovDataMap;
	// CVIPMapData<float> eventFovWeightDataMap;

	eventFovDataMap.Set( m_eventMapSize, "eventToFovMap_", fieldOfView.GetNumberOfVoxels() );

	CVIPMapData<unsigned int> eventFovDataMapNew;
	eventFovDataMapNew.Set( m_eventMapSize, "eventToFovMapNew_", fieldOfView.GetNumberOfVoxels() );


	/*
	if (CUserParameters::Instance()->GetUseSiddon() )
	{
		io_eventFovWeightDataMap.Set( CUserParameters::Instance()->GetEventMapSetSize()
				, "eventToFovWeightMap_"
				, fieldOfView.GetNumberOfVoxels() );
	}
	*/

	// 1. Read all the eventToFovMap files...
	int totalNumEvents = 0;
	// int numMapFiles = CUserParameters::Instance()->GetNumMapFiles();

	if (m_makeFakeFovMap)
	{
		MakeFakeFovMap( fieldOfView, eventFovDataMap );
		cout << "made a fake eventFovDataMap of size: " << eventFovDataMap.GetCurrentSize() << endl;
		totalNumEvents = eventFovDataMap.GetCurrentSize();
	}
	else
	{
		for (int imap = 0; imap < m_numMapFiles; imap++)
		{
			std::vector<unsigned int> fovidxVec;
			eventFovDataMap.Get(totalNumEvents, fovidxVec);
			totalNumEvents += eventFovDataMap.GetCurrentSize();
			/*
			if (CUserParameters::Instance()->GetUseSiddon() )
			{
				std::vector<float> fovweightVec;
				io_eventFovWeightDataMap.Get(totalNumEvents, fovweightVec);
			}
			*/
		}
		cout << "Total number of events (read from map files): " << totalNumEvents << endl;
	}

	int ievent = 0;
	C3Vector radiusVector(0.0, 0.0, 0.0);
	std::vector<unsigned int> fovidxVec;
	while (ievent < totalNumEvents)
	{
		if (m_makeFakeFovMap)
		{
			cout << endl;
			cout << "IEVENT: " << ievent << endl;
		}
		// Get the list of FOV indices for this event
		eventFovDataMap.Get(ievent, fovidxVec);
		std::vector<unsigned int> fovidxVecNew;

		PLANE_AXIS planeAxis = PLANE_AXIS_NONE;
		unsigned int firstidx = fovidxVec[0];
		C3Vector firstCentre; 
		fieldOfView.GetVoxelCentre(firstidx, firstCentre);
		unsigned int lastidx = fovidxVec[fovidxVec.size()-1];
		C3Vector lastCentre; 
		fieldOfView.GetVoxelCentre(lastidx, lastCentre);
		
		// We navigate along the direction in which the LOR has minimum gradient...
		GetLorPlaneAxis( firstidx, lastidx, fieldOfView, planeAxis );
		int ix, iy, iz;
		int prev_planeidx(-1);			

		// Loop over the FOV indices for this event
		for (int ij = 0; ij < fovidxVec.size(); ij++)
		{
			unsigned int fovidx = fovidxVec[ij];
			fieldOfView.GetVoxelIndices(fovidx, ix, iy, iz);
			
			C3Vector voxelCentre; 
			fieldOfView.GetVoxelCentre(fovidx, voxelCentre);

			// Add this voxel to the new list
			fovidxVecNew.push_back( fovidx );

			// 2. Check we are in the next plane
			int planeidx = (planeAxis == PLANE_AXIS_Z) ? iz 
					 	 : (planeAxis == PLANE_AXIS_Y) ? iy : ix;

			if (m_makeFakeFovMap)
			{
				cout << "fovidx: " << fovidx << ": [" << ix << "," << iy << "," << iz << "], with planeAxis: " << planeAxis << endl;
			}

			if (prev_planeidx == -1 || planeidx != prev_planeidx)
			{
				// check
				if (   !m_isIntersection			// intersection-sample might have two intersection points 
													// located at rather distant parts of the LOR
					&& prev_planeidx != -1 
					&& fabs(planeidx - prev_planeidx) > 1)
				{
					// Check that the current plane is right next to the previous plane
					//
					cout << planeidx << " " << prev_planeidx << endl;
					assert( false );
				}

				// Make tube radius
				if (planeAxis == PLANE_AXIS_Z)
					radiusVector.Set( m_tubeRadius, m_tubeRadius, 0.0);
				else if (planeAxis == PLANE_AXIS_Y)
					radiusVector.Set( m_tubeRadius, 0.0, m_tubeRadius);
				else
					radiusVector.Set( 0.0, m_tubeRadius, m_tubeRadius);

				// 3. Get coordinates of min and max corner
				C3Vector minPos( voxelCentre - radiusVector );		// left-bottom corner
// cout << "voxelCentre: " << voxelCentre << " radiusVector: " << radiusVector << " minPos: " << minPos << endl;
				if ( !fieldOfView.IsInBounds(minPos) )
				{
					if (minPos.GetX() < fieldOfView.GetLowerEdge().GetX())
						minPos.SetX( fieldOfView.GetLowerEdge().GetX() );
					if (minPos.GetY() < fieldOfView.GetLowerEdge().GetY())
						minPos.SetY( fieldOfView.GetLowerEdge().GetY() );
					if (minPos.GetZ() < fieldOfView.GetLowerEdge().GetZ())
						minPos.SetZ( fieldOfView.GetLowerEdge().GetZ() );
				}
// cout << "minPos revisited: " << minPos << endl;
				C3Vector maxPos( voxelCentre + radiusVector );		// right-top corner
				if ( !fieldOfView.IsInBounds(maxPos) )
				{
					if (maxPos.GetX() > fieldOfView.GetUpperEdge().GetX())
						maxPos.SetX( fieldOfView.GetUpperEdge().GetX() );
					if (maxPos.GetY() > fieldOfView.GetUpperEdge().GetY())
						maxPos.SetY( fieldOfView.GetUpperEdge().GetY() );
					if (maxPos.GetZ() > fieldOfView.GetUpperEdge().GetZ())
						maxPos.SetZ( fieldOfView.GetUpperEdge().GetZ() );
				}
					// TODO: actually, we should take the direction with an angle with respect to the radiusVector
					// i.e. the angle between the LOR and the axis perpendicular to the "z-plane"
					// This line - from the Lor to the tube edge, along the line perpendicular to z - is larger than radius.
					// Hence: If we take the "radiusVector", the effective tube is smaller than it should be.

				int imin1, imin2, imax1, imax2, idummin, idummax;
				int iminpx, iminpy, iminpz;
				int imaxpx, imaxpy, imaxpz;
				fieldOfView.GetVoxelIndices( minPos, iminpx, iminpy, iminpz );
				fieldOfView.GetVoxelIndices( maxPos, imaxpx, imaxpy, imaxpz );
				if (planeAxis == PLANE_AXIS_Z)
				{
					imin1 = iminpx; imin2 = iminpy; idummin = iminpz;
					imax1 = imaxpx; imax2 = imaxpy; idummax = imaxpz;
				}
				else if (planeAxis == PLANE_AXIS_Y)
				{
					imin1 = iminpx; imin2 = iminpz; idummin = iminpy;
					imax1 = imaxpx; imax2 = imaxpz; idummax = imaxpy;
				}
				else if (planeAxis == PLANE_AXIS_X)
				{
					imin1 = iminpy; imin2 = iminpz; idummin = iminpx;
					imax1 = imaxpy; imax2 = imaxpz; idummax = imaxpx;
				}
				assert (idummin == planeidx);
				assert (idummax == planeidx);

				if (m_makeFakeFovMap)
				{
					cout << "Additional: " ; // << endl;
					/*
					cout << "from minpos: " << minPos << " [" << iminpx << "," << iminpy << "," << iminpz << "]" 
						 << " and maxpos: " << maxPos << " [" << imaxpx << "," << imaxpy << "," << imaxpz << "]" 
						 << endl;
					*/
				}
				for (int i1 = imin1; i1 <= imax1; i1++)
				{
					for (int i2 = imin2; i2 <= imax2; i2++)
					{
						int index = (planeAxis == PLANE_AXIS_Z) ? fieldOfView.GetVoxelIndex( i1, i2, planeidx )
									: (planeAxis == PLANE_AXIS_Y) ? fieldOfView.GetVoxelIndex( i1, planeidx, i2 )
									: fieldOfView.GetVoxelIndex( planeidx, i1, i2 );
						// make sure we didn't add this voxel before...
/*
cout << "found an index: "  << index << endl; 
int kux, kuy, kuz;
fieldOfView.GetVoxelIndices(index, kux, kuy, kuz);
cout << index << ": [" << kux << "," << kuy << "," << kuz << "]  " << endl;
*/
						if ( !m_isIntersection
							&& index != fovidx 
							&& std::find(fovidxVecNew.begin(), fovidxVecNew.end(), index) != fovidxVecNew.end() )
						{
							// Check that index is not already added before.
							cout << "current index: " << index << " from fovidx: " << fovidx 
								 << " plane type: " << planeAxis << " plane idx: " << planeidx 
								 << endl;
							cout << "already added: " << endl;
							for (int ju = 0; ju < fovidxVecNew.size(); ju++)
								cout <<  "[" << ju << "]: " << fovidxVecNew[ju] << " ";
							cout << endl;
							assert( false );
						}
						else if (   index != fovidx 
								 && std::find(fovidxVecNew.begin(), fovidxVecNew.end(), index) == fovidxVecNew.end())
						{
							fovidxVecNew.push_back( index );
							if (m_makeFakeFovMap)
							{
								int kx, ky, kz;
								fieldOfView.GetVoxelIndices(index, kx, ky, kz);
								cout << index << ": [" << kx << "," << ky << "," << kz << "]  ";
							}

						}
					}	// loop over all additional voxels in direction 2
				} // loop over all additional voxels in direction 1
				if (m_makeFakeFovMap)
				{
					cout << " ENDadditional" << endl;
				}

				prev_planeidx = planeidx;
			}	// deal with new plane
		}	// loop over indices in Fov list of current event

		// Add this to the new map...
		eventFovDataMapNew.Add(ievent, fovidxVecNew);

		ievent++;
	} // loop over events

	eventFovDataMapNew.Finalize();
} 

// =========================================================================

void
GetLorPlaneAxis( int in_firstidx, int in_lastidx, const CVIPFieldOfView& in_fieldOfView, PLANE_AXIS& out_planeAxis )
{
	C3Vector voxelCentre1;
	in_fieldOfView.GetVoxelCentre(in_firstidx, voxelCentre1);
	C3Vector voxelCentre2;
	in_fieldOfView.GetVoxelCentre(in_lastidx, voxelCentre2);

	double delx = voxelCentre2.GetX() - voxelCentre1.GetX();
	double dely = voxelCentre2.GetY() - voxelCentre1.GetY();
	double delz = voxelCentre2.GetZ() - voxelCentre1.GetZ();

	if (fabs(delx) > fabs(dely) && fabs(delx) > fabs(delz))
		out_planeAxis = PLANE_AXIS_X;
	else if (fabs(dely) > fabs(delx) && fabs(dely) > fabs(delz))
		out_planeAxis = PLANE_AXIS_Y;
	else
		out_planeAxis = PLANE_AXIS_Z;
}
		
// =========================================================================

void
ReadParameters()
{
	std::string conffilename = "tor_parameters.conf";
	ifstream conffile(conffilename.c_str(), ios::in);
	if (!conffile.is_open())
	{
		cout << "Configure file not open: " << conffilename << endl;
		exit(1);
	}

	std::string tmpStr, tmpDum;
	double tmpVal;

	while (!conffile.eof())
    {
		conffile >> tmpStr >> tmpVal >> tmpDum;
		if (conffile.fail() && !conffile.eof())
		{
			cout << "ERROR! FORMAT ERROR in parameters file" << endl;
			exit(1);
		}
		if (!conffile.eof())
		{
			if (tmpStr == "eventMapSetSize")
			{
				m_eventMapSize = (int) tmpVal;
			}
			else if (tmpStr == "NumberOfMapFiles")
			{
				m_numMapFiles = (int) tmpVal;
			}
			else if (tmpStr == "TubeRadius")
			{
				m_tubeRadius = tmpVal;
			}

			// ---------------------------------------------------------
			else if (!conffile.eof())
			{
				cout << "ERROR! Unknown entry in configure file: " << tmpStr << endl;
				exit(1);
			}
		}
	}
}

// =========================================================================

void
MakeFakeFovMap( const CVIPFieldOfView& in_fieldOfView, CVIPMapData<unsigned int>& io_eventFovDataMap )
{
	int nbinsX = in_fieldOfView.GetNumVoxelsX();
	int nbinsY = in_fieldOfView.GetNumVoxelsY();
	int nbinsZ = in_fieldOfView.GetNumVoxelsZ();

	// ievent 1 (= LOR 1)
	{
		int iy = (int) (0.5 * nbinsY);
		int iz = 0;
	
		std::vector<unsigned int> fovidxVec;
		for (int ix = 0; ix < nbinsX; ix++)
		{
			int index = in_fieldOfView.GetVoxelIndex( ix, iy, iz );
			fovidxVec.push_back( index );
		}
	
		io_eventFovDataMap.Add(0, fovidxVec);
	}
	// ievent 2 (= LOR 2)
	{
		int ix = (int) (0.5 * nbinsX);
		int iz = 0;
	
		std::vector<unsigned int> fovidxVec;
		for (int iy = 0; iy < nbinsY; iy++)
		{
			int index = in_fieldOfView.GetVoxelIndex( ix, iy, iz );
			fovidxVec.push_back( index );
		}
	
		io_eventFovDataMap.Add(1, fovidxVec);
	}
	// ievent 3 (= LOR 3)
	{
		int ix = 2;
		int iz = 0;
	
		std::vector<unsigned int> fovidxVec;
		for (int iy = 0; iy < nbinsY; iy++)
		{
			int index = in_fieldOfView.GetVoxelIndex( ix, iy, iz );
			fovidxVec.push_back( index );
			if (iy%4 == 0)
				ix++;
		}
	
		io_eventFovDataMap.Add(2, fovidxVec);
	}

	io_eventFovDataMap.Finalize();
}









