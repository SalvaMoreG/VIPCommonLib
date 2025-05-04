#include "CVIPImage.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

/*
Programname: ConvertImgToWeighted3DTree
This code reads a binary LM-OSEM output file (e.g. "out_binary_amide_ITER4.img"), input:
    - out_binary_amide_ITER4.img
    - osem_fov_parameters.conf
and converts the information into output:
    - a root file with the following contents, for each 3D image bin:
        x, y, z, weight (i.e. pixel-content)
    - a root file with histograms: 'histograms.root'
*/

int main()
{
    cout << "Programname: ConvertImgToWeighted3DTree" << endl;
    cout << "This code reads a binary LM-OSEM output file, input: " << endl;
    cout << "   - out_binary_amide_ITER4.img" << endl;
    cout << "   - osem_fov_parameters.conf" << endl;
    cout << "and converts the information into output:" << endl;
    cout << "   - a root file with the following contents, for each 3D image bin: " << endl;
    cout << "       x, y, z, weight (i.e. pixel-content)" << endl;
    cout << "   - a root file with histograms: 'histograms.root'" << endl;
    cout << endl;

    // FOV
    string filename;
    CVIPFieldOfView fieldOfView;
    cout << "Give FOV parameters filename: " << endl;
    cin >> filename;
    fieldOfView.Initialize( filename.c_str() );
    int nVoxels = fieldOfView.GetNumberOfVoxels();
    cout << "nVoxels: " << nVoxels << endl;

    // Create Image
    cout << "Give binary image filename: " << endl;
    cin >> filename;

    // Reading
    CVipImage theImage;
    theImage.Read( filename, fieldOfView, FFORMAT_BINARY, DFORMAT_FLOAT);

    // And writing a check image (should be same as input)
    std::string out_filename( "check_image.img" );
    theImage.Write( out_filename, fieldOfView, FFORMAT_BINARY );

    // checking histograms
    double maxX = fieldOfView.GetSize().GetX()/2.0;
    TH1D* hX = new TH1D("hX", "All along X", fieldOfView.GetNumVoxelsX(), -1.0*maxX, maxX);
    double maxY = fieldOfView.GetSize().GetY()/2.0;
    TH1D* hY = new TH1D("hY", "All along Y", fieldOfView.GetNumVoxelsY(), -1.0*maxY, maxY);
    double maxZ = fieldOfView.GetSize().GetZ()/2.0;
    TH1D* hZ = new TH1D("hZ", "All along Z", fieldOfView.GetNumVoxelsZ(), -1.0*maxZ, maxZ);

    // the 3D TTree:
    TTree* the3DTree = new TTree("the3DTree", "3D");
    double tt_x, tt_y, tt_z, tt_value;
    if (the3DTree)
    {
        the3DTree->Branch("x", &tt_x);
        the3DTree->Branch("y", &tt_y);
        the3DTree->Branch("z", &tt_z);
        the3DTree->Branch("WEIGHT", &tt_value);
    }

    C3Vector coordinates;
    float checksum(0.0);
    for ( int ivox = 0; ivox < nVoxels; ivox++ )
    {
		tt_value = 0.0;
        if ( fieldOfView.GetVoxelCentre(ivox, coordinates) )
        {
        	tt_value = theImage.GetVoxelValue( ivox );
            // ASCII output
            /*
            out_file << std::left << std::setw(12) << coordinates.GetX()
                        << std::left << std::setw(12) << coordinates.GetY()
                        << std::left << std::setw(12) << coordinates.GetZ()
                        << value << endl;
            */

            hX->Fill(coordinates.GetX(), tt_value);
            hY->Fill(coordinates.GetY(), tt_value);
            hZ->Fill(coordinates.GetZ(), tt_value);

            tt_x = coordinates.GetX();
            tt_y = coordinates.GetY();
            tt_z = coordinates.GetZ();
            the3DTree->Fill();
        }
    }

    TFile* myfile = TFile::Open("histograms.root","RECREATE");
    hX->Write();
    hY->Write();
    hZ->Write();
    myfile->Close();

    TFile aFile("3D_TREE.root_NEW", "RECREATE");
    the3DTree->Write();

}

