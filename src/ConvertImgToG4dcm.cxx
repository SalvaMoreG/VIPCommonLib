 
#include "CVIPImage.h"

#include "TFile.h"
#include "TH1D.h"
 
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

/*
Programname: ConvertImgToG4dcm
This code reads a binary LM-OSEM output file (e.g. "out_binary_amide_ITER4.img"), input: 
	- out_binary_amide_ITER4.img
	- osem_fov_parameters.conf
and converts the information into output: 
	- a "g4dcm" file: check_image_ascii.g4dcm
	- a root file with projected line-profiles: line_profile_histograms.root
	- a check image file: check_image.img
*/

int main()
{
	cout << "Programname: ConvertImgToG4dcm" << endl;
	cout << "This code reads a binary LM-OSEM output file, input: " << endl;
	cout << "	- out_binary_amide_ITER4.img" << endl;
	cout << "	- osem_fov_parameters.conf" << endl;
	cout << "and converts the information into output:" << endl; 
	cout << "	- a 'g4dcm' file: check_image_ascii.g4dcm" << endl;
	cout << "	- a root file with projected line-profiles: line_profile_histograms.root" << endl;
	cout << "	- a check image file: check_image.img" << endl;
	cout << endl;
 
    cout << "Give image filename: " << endl;
    string filename; cin >> filename;
    
	CVIPFieldOfView fieldOfView;
	fieldOfView.Initialize( "osem_fov_parameters.conf" );    

	// Create Image
	int nVoxels = fieldOfView.GetNumberOfVoxels();
    cout << "nVoxels: " << nVoxels << endl;

	cout << "New g4dcm file will have #slices (in Z): " << fieldOfView.GetNumVoxelsZ() << endl;
	cout << "i.e., we will loop from 0 to " << fieldOfView.GetNumVoxelsZ()-1 << endl;
	int numVoxelsZ( fieldOfView.GetNumVoxelsZ() );
	int minZbin(0), maxZbin(numVoxelsZ-1);
	double binsizeZ = (fieldOfView.GetUpperEdge().GetZ() - fieldOfView.GetLowerEdge().GetZ());
	binsizeZ = binsizeZ / fieldOfView.GetNumVoxelsZ(); 
	cout << "bin size (in Z): " << binsizeZ << endl;
	cout << "Give min bin and max bin index (in Z) (0 " << maxZbin << " means no change)" << endl;
	cin >> minZbin >> maxZbin;
	if (minZbin != 0 && maxZbin != numVoxelsZ-1)
		numVoxelsZ = (maxZbin - minZbin) + 1;

    cout << "New g4dcm file will have #slices (in Z): " << numVoxelsZ << endl;
    cout << "i.e., we will loop from " << minZbin << " to " << maxZbin << endl;

    // Reading
    CVipImage theImage;	
    theImage.Read( filename, fieldOfView, FFORMAT_BINARY, DFORMAT_FLOAT);
    
    // And writing a check image (should be same as input)
    std::string out_filename( "check_image.img" );
    theImage.Write( out_filename, fieldOfView, FFORMAT_BINARY ); 
    
    // std::string out_filename_ascii( "check_image_ascii.txt" );
    // ofstream out_file( out_filename_ascii.c_str() );
    
    std::string out_filename_g4dcm( "check_image_ascii.g4dcm" );
    ofstream out_file_g4dcm( out_filename_g4dcm.c_str() );
    
    double maxX = fieldOfView.GetSize().GetX()/2.0; 
    TH1D* hX = new TH1D("hX", "All along X", fieldOfView.GetNumVoxelsX(), -1.0*maxX, maxX);         
    double maxY = fieldOfView.GetSize().GetY()/2.0; 
    TH1D* hY = new TH1D("hY", "All along Y", fieldOfView.GetNumVoxelsY(), -1.0*maxY, maxY); 
    double maxZ = fieldOfView.GetSize().GetZ()/2.0; 
    TH1D* hZ = new TH1D("hZ", "All along Z", fieldOfView.GetNumVoxelsZ(), -1.0*maxZ, maxZ);

    C3Vector coordinates;
    float checksum(0.0);
    for ( int ivox = 0; ivox < nVoxels; ivox++ )
    {
        float value = theImage.GetVoxelValue( ivox );
        
        if ( fieldOfView.GetVoxelCentre(ivox, coordinates) )
        {
            // ASCII output
            /*
            out_file << std::left << std::setw(12) << coordinates.GetX()
                        << std::left << std::setw(12) << coordinates.GetY() 
                        << std::left << std::setw(12) << coordinates.GetZ() 
                        << value << endl;
            */
            hX->Fill(coordinates.GetX(), value);
            hY->Fill(coordinates.GetY(), value);
            hZ->Fill(coordinates.GetZ(), value);

            if (value > 0.0 && coordinates.GetZ() == 0.0 )
            {
                /*
                cout << std::left << std::setw(8) << coordinates.GetX()
                        << std::left << std::setw(8) << coordinates.GetY() 
                        << std::left << std::setw(8) << coordinates.GetZ() 
                        << std::left << std::setw(8) << value << endl;  
                */
                checksum += value;
            }
        }
    }
    // cout << "checksum: " << checksum << endl;

    int ix, iy, iz, ivox;  
    double value;
    out_file_g4dcm << fieldOfView.GetNumVoxelsX() << " "
                   << fieldOfView.GetNumVoxelsY() << " "
                   // << fieldOfView.GetNumVoxelsZ() << endl;
				   << numVoxelsZ << endl;
    out_file_g4dcm << fieldOfView.GetLowerEdge().GetX() << " " 
                   << fieldOfView.GetUpperEdge().GetX() << endl;
    out_file_g4dcm << fieldOfView.GetLowerEdge().GetY() << " " 
                   << fieldOfView.GetUpperEdge().GetY() << endl;    
	double lowerZ = fieldOfView.GetLowerEdge().GetZ() + minZbin * binsizeZ;
	double upperZ = lowerZ + numVoxelsZ * binsizeZ;
    out_file_g4dcm << lowerZ << " " << upperZ << endl;

	binsizeZ = binsizeZ / fieldOfView.GetNumVoxelsZ(); 
                   
    // for (iz = 0; iz < fieldOfView.GetNumVoxelsZ(); iz++)
    for (iz = minZbin; iz < maxZbin+1; iz++)
    {
        for (iy = 0; iy < fieldOfView.GetNumVoxelsY(); iy++)
        {
            for (ix = 0; ix < fieldOfView.GetNumVoxelsX(); ix++)
            {
                ivox = fieldOfView.GetVoxelIndex( ix, iy, iz );
                value = theImage.GetVoxelValue( ivox );
                out_file_g4dcm << value << " ";
            } 
            out_file_g4dcm << endl;
        }
    }
    
    TFile* myfile = TFile::Open("line_profile_histograms.root","RECREATE");
    hX->Write();
    hY->Write();
    hZ->Write();
    myfile->Close();
}

