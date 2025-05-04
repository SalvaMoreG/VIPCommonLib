
#include "CVIPMatrix.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
	cout << "Give filename" << endl;
	string fname; cin >> fname;
	ifstream infile( fname.c_str() );
	if ( !infile.is_open() )
	{
		cout << "File not found: " << fname << endl; 
		return -1;
	}
	
	cout << "isCT (0 = no/1 = yes)?" << endl;
	bool isCT; cin >> isCT;
	ofstream outfile;
	if (isCT)
	{
		cout << "CT Image" << endl;
		outfile.open("newCT_file.g4dcm_FLIPPED_NEW");
	} else {
		cout << "PET Image" << endl;
		outfile.open("newPET_file.g4dcm_FLIPPED_NEW");
	}
	
    if (isCT)
    {
        string str;
        for (int i = 0; i < 5; i++)
        {
            getline(infile, str);
            outfile << str << endl;
        }
    }
	int nX, nY, nZ;
	infile >> nX >> nY >> nZ;
	double minX, maxX, minY, maxY, minZ, maxZ;
	infile >> minX >> maxX;
	infile >> minY >> maxY;
	infile >> minZ >> maxZ;
    
    outfile << nX << " " << nY << " " << nZ << endl;
    outfile << minX << " " << maxX << endl;
    outfile << minY << " " << maxY << endl;
    outfile << minZ << " " << maxZ << endl;
            
	double posX, posY, posZ;
    CVIPMatrix* matrix = new CVIPMatrix ( nY, nX );
	if (isCT)
	{
		int matidx;
		for (int iZ = 0; iZ < nZ; iZ++)
		{
            // Read all XY values
			for (int iY = 0; iY < nY; iY++)
			{
				for (int iX = 0; iX < nX; iX++)
                {
					infile >> matidx;
                    matrix->SetElement( iY, iX, matidx );
                }
			}
			// Write all XY values, but upside down
			for (int iY = nY-1; iY >= 0; iY--)
			{
				for (int iX = 0; iX < nX; iX++)
                {
                    double tmp = matrix->GetElement( iY, iX );
                    outfile << (int) tmp << " ";
                }
                outfile << endl;
			}
		}
	}        
	
	// 
	
	double value; // CT: matdensity or PET: activity
	for (int iZ = 0; iZ < nZ; iZ++)
	{
		for (int iY = 0; iY < nY; iY++)
		{
			for (int iX = 0; iX < nX; iX++)
			{
				infile >> value;
                matrix->SetElement( iY, iX, value );
            }
        }
        // Write all XY values, but upside down
        for (int iY = nY-1; iY >= 0; iY--)
        {
            for (int iX = 0; iX < nX; iX++)
            {
                value = matrix->GetElement( iY, iX );
                outfile << value << " ";
                if (isCT)
                {
                    if ( (iX+1)%8 == 0) outfile << endl;
                }                
            }
            if (!isCT) outfile << endl;
        }        
    }

	return 0;
}

