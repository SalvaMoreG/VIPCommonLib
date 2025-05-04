 
#include "CVIPFieldOfView.h"
#include "CVIPImage.h"
 
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

void
WriteImageQualityFile( const CVIPFieldOfView& in_fieldOfView );

int main( int argc, char *argv[] )
{
	int ichoice;
	cout << "Make your choice: " << endl;
	cout << "(1) From FOV parameters file to Img Quality file 'imgqual_parameters.conf_NEW'" << endl;
	cout << "(2) From Img Quality file to FOV parameters 'fov_parameters.conf_NEW'" << endl;
	cin >> ichoice;

	if (ichoice < 1 || ichoice > 2) 
	{
		cout << "you made the wrong choice: " << ichoice << endl;
		return -1;
	}
	
	if (ichoice == 1)
	{
		CVIPFieldOfView fieldOfView;
		string parfname("osem_fov_parameters.conf");
		cout << "FOV parameters input filename (or type 0 if using default: " << parfname << " )" << endl;
		string tmpStr;
		cin >> tmpStr;
		if (tmpStr.size() > 1)
		{
			parfname = tmpStr;
		}
		cout << "Opening file: " << parfname << endl;
		fieldOfView.Initialize( parfname );

		WriteImageQualityFile( fieldOfView );
	}
}

void
WriteImageQualityFile( const CVIPFieldOfView& in_fieldOfView )
{
/*
FovBinsX                        161          //
FovBinsY                        161          //
FovBinsZ                        1         //

FovMinX                        -20.125        //
FovMinY                        -20.125        //
FovMinZ                        -12.5      //

FovVoxelSizeX                   0.25         //
FovVoxelSizeY                   0.25         //
FovVoxelSizeZ                   25.0         //
*/
	ofstream outfile("imgqual_parameters.conf_NEW");

	outfile << "FovBinsX                        " << in_fieldOfView.GetNumVoxelsX() << "         //" << endl;
	outfile << "FovBinsY                        " << in_fieldOfView.GetNumVoxelsY() << "         //" << endl;
	outfile << "FovBinsZ                        " << in_fieldOfView.GetNumVoxelsZ() << "         //" << endl;
	outfile << endl;
	outfile << "FovMinX                         " << in_fieldOfView.GetLowerEdge().GetX() << "         //" << endl;
	outfile << "FovMinY                         " << in_fieldOfView.GetLowerEdge().GetY() << "         //" << endl;
	outfile << "FovMinZ                         " << in_fieldOfView.GetLowerEdge().GetZ() << "         //" << endl;
	outfile << endl;
	outfile << "FovVoxelSizeX                   " << in_fieldOfView.GetVoxelSize().GetX() << "         //" << endl;
	outfile << "FovVoxelSizeY                   " << in_fieldOfView.GetVoxelSize().GetY() << "         //" << endl;
	outfile << "FovVoxelSizeZ                   " << in_fieldOfView.GetVoxelSize().GetZ() << "         //" << endl;
}



