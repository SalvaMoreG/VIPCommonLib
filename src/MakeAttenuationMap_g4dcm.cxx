#include "CVIPImage.h"
#include "CVIPFieldOfView.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

using namespace std;

void parseAttenuationConfFile( std::vector<double>& out_attcoeff_pergram );
bool is_even( int in_value );

int main( int argc, char* argv[] )
{
	// PARSING
    int m_zslices(0);
    double m_zsize(0);
	if (argc > 1)
    for (int iarg = 1; iarg < argc; iarg++)
	{
		string flag = argv[iarg];
		if ( flag == "-paddingZ" )
		{
            if ( iarg+1 >= argc )
            {
                    cout << "Wrong number of arguments. Usage: " << endl;
                    cout << argv[0] << " -paddingZ <numberOfTotalSlices>" << endl;
                    exit(1);
            }
            iarg++;
			m_zslices = atoi ( argv[iarg] );
		}
		else if ( flag == "-sizeZ" )
		{
            if ( iarg+1 >= argc )
            {
                    cout << "Wrong number of arguments. Usage: " << endl;
                    cout << argv[0] << " -sizeZ <sizeZ>" << endl;
                    exit(1);
            }
            iarg++;
			m_zsize = atof ( argv[iarg] );
		}
		else
		{
			cout << "Allowed flags are: " << endl;
			cout << argv[0] << " -paddingZ <numberOfTotalSlices>" << endl;
            cout << argv[0] << " -sizeZ <sizeZ>" << endl;
			return 1;
		}
	}

	cout << "give g4dcm filename: " << endl;
	string filename; cin >> filename;
	ifstream infile(filename.c_str());

	std::vector<double> attcoeff_pergram;
	parseAttenuationConfFile( attcoeff_pergram );

	int numMaterials(0);
	infile >> numMaterials;
	cout << "number of materials: " << numMaterials << endl;
	int index;
	string nameStr;
	for (int idx = 0; idx < numMaterials; idx++)
	{
		infile >> index >> nameStr;
		cout << "[" << idx << "]: " << nameStr << ": " << attcoeff_pergram[idx] << endl;
		if (index != idx)
			cout << "ERROR! index: " << index << " idx: " << idx << endl;
		assert( index == idx );
	}

	// Read CT image size
	int nbinsX, nbinsY, nbinsZ;
	infile >> nbinsX >> nbinsY >> nbinsZ;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	infile >> xmin >> xmax;
	infile >> ymin >> ymax;
	infile >> zmin >> zmax;
	// /*
	cout << "X: " << nbinsX << "[" << xmin << ";" << xmax << "]" << endl;
	cout << "Y: " << nbinsY << "[" << ymin << ";" << ymax << "]" << endl;
	cout << "Z: " << nbinsZ << "[" << zmin << ";" << zmax << "]" << endl;
	// */

    int half_padding_bins_Z(0);
    if (nbinsZ > m_zslices)
    {
        cout << "CT image has more slices than user-selected Z padding. Using CT image slices: "
             << nbinsZ << endl;
        m_zslices = nbinsZ;
    }
    else if (nbinsZ < m_zslices)
    {
        if (is_even(nbinsZ) != is_even(m_zslices))
        {
            cout << "CT image slices and user-selected Z padding should both be either even or odd"
                 << " CT image slices: " << nbinsZ << " user Z padding: " << m_zslices << endl;
            cout << "Using CT image slices: " << nbinsZ << endl;
            m_zslices = nbinsZ;
        }
        else
        {
            int padding_bins_Z = m_zslices - nbinsZ;
            if ( !is_even(padding_bins_Z) )
            {
                // should never happen
                cout << "ERROR! Additional padding bins is not an even number: "
                     << padding_bins_Z << endl;
                return 1;

            }
            double binsizeZ = (zmax - zmin)/nbinsZ;
            if (m_zsize > 0.0)
                binsizeZ = m_zsize;
            double midZ = zmin + (zmax - zmin)/2.0;
            zmin = midZ - 0.5 * m_zslices * binsizeZ;
            zmax = midZ + 0.5 * m_zslices * binsizeZ;

            half_padding_bins_Z = padding_bins_Z/2;
        }
    }

	CVIPFieldOfView fieldOfView;
    fieldOfView.SetFieldOfView( nbinsX, xmin, xmax, nbinsY, ymin, ymax, m_zslices, zmin, zmax );
	int ntotbins = fieldOfView.GetNumberOfVoxels();

	CVipImage materialImage( ntotbins );
	int materialID;
    // for (int iz = 0; iz < nbinsZ; iz++) {
	for (int iz = half_padding_bins_Z; iz <  half_padding_bins_Z + nbinsZ; iz++) {
		for (int iy = 0; iy < nbinsY; iy++) {
			for (int ix = 0; ix < nbinsX; ix++)
			{
				infile >> materialID;
				index =	fieldOfView.GetVoxelIndex( ix, iy, iz );
				// cout << "[" << ix << "," << iy << "," << iz << "] index: " << index
                //     << " matID: " << materialID << endl;
				materialImage.SetVoxelValue( index, materialID );
			}
		}
	}

	CVipImage attenuationImage( ntotbins );
	double density, attenuation;
	const double percm2permm(0.1); // in OSEM we will work with attenuation per mm
	// for (int iz = 0; iz < nbinsZ; iz++) {
    for (int iz = half_padding_bins_Z; iz <  half_padding_bins_Z + nbinsZ; iz++)
    {
		for (int iy = 0; iy < nbinsY; iy++)
		{
			for (int ix = 0; ix < nbinsX; ix++)
			{
				infile >> density;
				index =	fieldOfView.GetVoxelIndex( ix, iy, iz );
				materialID = (int) materialImage.GetVoxelValue( index );
				attenuation = density * attcoeff_pergram[materialID] * percm2permm;
				// cout << "[" << ix << "," << iy << "," << iz << "] index: " << index
                //     << " matID: " << materialID << " att: " << attenuation << endl;
				attenuationImage.SetVoxelValue( index, attenuation );
			}
		}
	}

	string out_fname("out_attmap.img_NEW");
	attenuationImage.Write( out_fname, fieldOfView, FFORMAT_BINARY );

	cout << "Attenuation image map ready" << endl;


    cout << endl;
    cout << "att_fov_parameters.conf: " << endl;
    cout << "minX           " << xmin << "    // " << endl;
    cout << "maxX           " << xmax << "    // " << endl;
    cout << "nVoxelsX       " << nbinsX << "    // " << endl;
    cout << "binsize: " << (xmax - xmin)/nbinsX << endl;
    cout << "minY           " << ymin << "    // " << endl;
    cout << "maxY           " << ymax << "    // " << endl;
    cout << "nVoxelsY       " << nbinsY << "    // " << endl;
    cout << "binsize: " << (ymax - ymin)/nbinsY << endl;
    cout << "minZ           " << zmin << "    // " << endl;
    cout << "maxZ           " << zmax << "    // " << endl;
    cout << "nVoxelsZ       " << m_zslices << "    // " << endl;
    cout << "binsize: " << (zmax - zmin)/m_zslices << endl;

	return 0;
}

void parseAttenuationConfFile( std::vector<double>& out_attcoeff_pergram )
{
	// std::vector<double> attcoeff_pergram{ 0.0, 0.09687, 0.0964, 0.0902 };
	// out_attcoeff_pergram = attcoeff_pergram;

	string filename("attenuation_pergram.conf");
	ifstream infile( filename.c_str() );
	if ( !infile.is_open() )
	{
		cout << "File not open: " << filename << endl;
		exit(1);
	}
	int index;
	string nameStr;
	double attcoeff_pergram;
	while ( !infile.eof() )
	{
		infile >> index >> nameStr >> attcoeff_pergram;
		out_attcoeff_pergram.push_back(attcoeff_pergram);
	}
}

bool is_even( int in_value )
{
    return (in_value == 2 * (in_value/2));
}



