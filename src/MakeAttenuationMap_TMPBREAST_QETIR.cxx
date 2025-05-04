#include "CVIPFieldOfView.h"
#include "CVIPImage.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

const double mu_value_material = 0.00968;

int main( int argc, char *argv[] )
{
	CVIPFieldOfView fieldOfView;
	fieldOfView.Initialize( "fov_parameters.conf" );

	string fname;
	cout << "give starting image (*.img) filename" << endl;
	cin >> fname;

	CVipImage image;
	FILE_FORMAT fileformat( FFORMAT_BINARY );
	DATA_FORMAT dataformat( DFORMAT_FLOAT );
	image.Read( fname, fieldOfView, fileformat, dataformat);

	int nvoxels = image.GetSize();

	int nvoxels_check = fieldOfView.GetNumberOfVoxels();
	if ( nvoxels != nvoxels_check )
	{
		cout << "Image size: " << nvoxels << endl;
		cout << "FOV size: " << nvoxels_check << endl;
	}
	assert( nvoxels == nvoxels_check );

	cout << "nvoxels: " << nvoxels << endl;
	cout << "give minimum signal threshold" << endl;
	double threshold;
	cin >> threshold;

	int ix, iy, iz;
	int numX = fieldOfView.GetNumVoxelsX();
	int numY = fieldOfView.GetNumVoxelsY();
	for (int iv = 0; iv < nvoxels; iv++) 
	{
		fieldOfView.GetVoxelIndices( iv, ix, iy, iz);

		double val = image.GetVoxelValue( iv );
		if (val < threshold || ix < 2 || ix > numX-3 || iy < 2 || iy > numY-3)
			image.SetVoxelValue( iv, 0.0 );
		else
			image.SetVoxelValue( iv, mu_value_material );

		C3Vector coords;
		fieldOfView.GetVoxelCentre( iv, coords );
		double radiusXY2 = coords.GetX()*coords.GetX() + coords.GetY()*coords.GetY();
		double distZ = -1.0 * (coords.GetZ() - 20.0);
		double radius = sqrt(radiusXY2 + distZ*distZ);
		if (radius < 70.0 && coords.GetZ() < 20)
		{
			image.SetVoxelValue( iv, mu_value_material );
		}
	}

	string fname_new( "ATTENUATION_map.img_NEW" );
	image.Write( fname_new, fieldOfView, fileformat );

}



