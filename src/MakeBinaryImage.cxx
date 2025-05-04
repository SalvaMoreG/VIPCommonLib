
#include "CVIPImage.h"
#include "CVIPFieldOfView.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

int main( int argc, char *argv[] )
{
	CVIPFieldOfView fieldOfView;
    string infilename;
    cout << "give name of 'FOV' geometry file" << endl;
    cin >> infilename;
	fieldOfView.Initialize( infilename.c_str() );

	int nvoxels = fieldOfView.GetNumberOfVoxels();

    /*
	cout << "give standard cell contents" << endl;
	double value;
	cin >> value;
    */
    cout << "give input filename (make sure #entrees = " << nvoxels << ")" << endl;
    
    cin >> infilename;
    ifstream infile( infilename.c_str() );
    if ( !infile.is_open() )
    {
        cout << "File not found: " << infilename << endl;
        return -1;
    }

	int ix, iy, iz;
	int numX = fieldOfView.GetNumVoxelsX();
	int numY = fieldOfView.GetNumVoxelsY();
	cout << "nvoxels: " << nvoxels << " so, X x Y =" << numX << " x " << numY << endl;
	CVipImage image( nvoxels );
    
    int iv(0);
    double value(0.0);
	// for (int iv = 0; iv < nvoxels; iv++) 
    while ( !infile.eof() && iv < nvoxels )
	{
        infile >> value;
		image.SetVoxelValue( iv, value );

        /*
		fieldOfView.GetVoxelIndices( iv, ix, iy, iz);        
		if ( iv == 100 || (ix == (int) (0.5 * numX) && iy == (int) (0.5 * numY)) )
		{
cout << "got a special one, at: [" << ix << "," << iy << "," << iz << "]" << endl;
			image.SetVoxelValue( iv, 0.5*value );
		}
		*/
        iv++;
	}

	string fname_new( "BINARY_IMAGE.img_NEW" );
	FILE_FORMAT fileformat( FFORMAT_BINARY );
	image.Write( fname_new, fieldOfView, fileformat );

}
