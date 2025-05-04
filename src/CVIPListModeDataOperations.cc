
#include "CVIPListModeDataOperations.h"
#include "VIPconstants.h"

#include "CVIPRandom.h"

#include <iostream>

using namespace std;

// TODO: find a way to get rid of the following hard-coded stuff

// PET
const double m_voxelSizeY = 1.0; // mm TODO: of course, this should not be hardcoded!!!

// PEM
const double m_voxelSizeX_PEM = 2.0; // mm TODO: of course, this should not be hardcoded!!!
const double m_voxelSizeZ_PEM = 1.0; // mm TODO: of course, this should not be hardcoded!!!

void
CVIPListModeDataOperations::TangentialBlurring( C3Vector& io_position )
{
    // cout << "in the beginning: " << io_position << endl;
	// 1. Calculate angle for this event
	double angle(0);
	if (io_position.GetX() != 0 )
	{
		angle = atan( io_position.GetY()/io_position.GetX() );		// TODO: test for all quadrants!
	}
	else
	{
		angle = io_position.GetY() > 0 ? 0.5 * kPI : -0.5 * kPI;
	}

	// 2. Rotate to 0 degrees
	m_matrix.Set( Axis_Z, -1.0 * angle);
	io_position = m_matrix * io_position;

	// 3. Set random position along Y
	CVIPRandomUniform uniformDist(0.0, 1.0);
	double shot = m_voxelSizeY * uniformDist.GetNewValue();

	double newy = shot - 0.5 * m_voxelSizeY;

	io_position.SetY( newy );

	// 4. Rotate back to original position
	m_matrix.Set( Axis_Z, angle);
	io_position = m_matrix * io_position;
}

// ==============================================================

void
CVIPListModeDataOperations::TangentialBlurringPEM( C3Vector& io_position ) const
{
	// 1. Set random position along Z
	CVIPRandomUniform uniformDist(0.0, 1.0);
	double dZ = m_voxelSizeZ_PEM * uniformDist.GetNewValue();
	double dX = m_voxelSizeX_PEM * uniformDist.GetNewValue();

    // double newval = io_position.GetZ() + CLHEP::RandFlat::shoot( m_voxelSizeZ_PEM ) - 0.5 * m_voxelSizeZ_PEM;
	double newval = io_position.GetZ() + dZ - 0.5 * m_voxelSizeZ_PEM;
	io_position.SetZ( newval );

	// 2. Set random position along X
	// newval = io_position.GetX() + CLHEP::RandFlat::shoot( m_voxelSizeX_PEM ) - 0.5 * m_voxelSizeX_PEM;
    newval = io_position.GetX() + dX - 0.5 * m_voxelSizeX_PEM;
	io_position.SetX( newval );
}


/*
cout << "TESTING ANGLES <<<<<<" << endl;
double x, y;
for (int k = 0; k < 4; k++)
{
	if (k == 0)
	{
		x = 3.0; y = 1.0;
	}
	else if (k == 1)
	{
		x = -1.0; y = 3.0;
	}
	else if (k == 2)
	{
		x = -3.0; y = -1.0;
	}
	else if (k == 3)
	{
		x = 1.0; y = -3.0;
	}
	angle = atan( y/x );		// TODO: test for all quadrants!
cout << "		angle: " << angle*180/kPI << endl;
}
cout << " >>>>>>>>> " << endl;
cout << endl;
*/





