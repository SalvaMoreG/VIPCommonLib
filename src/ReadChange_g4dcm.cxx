
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

bool k_SetOutsideSelectedRegionToZero(false);

int main(int argc, char* argv[])
{
	if (argc > 1)
	{
		string flag = argv[1];
		if (flag == "-settozero")
		{
			k_SetOutsideSelectedRegionToZero = true;
			cout << "All pixels outside the selected region are set to 0" << endl;
		}
		else
		{
			cout << argv[0] << " -settozero" << endl;
			return 0;
		}
	}

	string fname;
	cout << "give filename" << endl;
	cin >> fname;
	ifstream ffile( fname.c_str() );
	if ( !ffile.is_open() )
	{
		cout << "file does not exist: " << fname << endl;
		return 1;
	}

	cout << "isCT (0 = no/1 = yes)?" << endl;
	bool isCT; cin >> isCT;
	ofstream outfile;
	if (isCT)
	{
		cout << "CT Image" << endl;
		outfile.open("newCT_file.g4dcm_NEW");
	} else {
		cout << "PET Image" << endl;
		outfile.open("newPET_file.g4dcm_NEW");
	}
	
	// Read header
	if (isCT)
	{
		int numMaterials(0);
		int dumInt;
		string dumStr;
		ffile >> numMaterials;
		outfile << numMaterials << endl;
		for (int i = 0; i < numMaterials; i++)
		{
			ffile >> dumInt >> dumStr;
			outfile << dumInt << " " << dumStr << endl;
		}
	}

	int nX, nY, nZ;
	ffile >> nX >> nY >> nZ;
	double minX, maxX, minY, maxY, minZ, maxZ;
	ffile >> minX >> maxX;
	ffile >> minY >> maxY;
	ffile >> minZ >> maxZ;
	double binSizeX = (maxX - minX)/nX;
	double binSizeY = (maxY - minY)/nY;
	double binSizeZ = (maxZ - minZ)/nZ;

	cout << "Give Radius of PET Scanner" << endl;	
	double radius; cin >> radius;
	double maxSize = radius/(std::sqrt(2));
	cout << "maxSize in X and Y = [ " << -maxSize << " : " << maxSize << " ] " << endl;

	double shiftY(0.0); 
	if ( !k_SetOutsideSelectedRegionToZero )
	{
		cout << "Give shift in Y?" << endl;
		cin >> shiftY;
	}

	int newNX = int(2.0*maxSize/binSizeX);
	int newNY = int(2.0*maxSize/binSizeY);
	int newNZ = nZ;
	double newMaxX = 0.5 * newNX * binSizeX;
	double newMinX = -1.0 * newMaxX;
	double newMaxY = 0.5 * newNY * binSizeY;
	double newMinY = -1.0 * newMaxY;
	newMaxY += shiftY;
	newMinY += shiftY;
	double newMinZ(minZ);
	double newMaxZ(maxZ);
	cout << newNX << " " << newNY << " " << newNZ << endl;	
	cout << newMinX << " " << newMaxX << endl;
	cout << newMinY << " " << newMaxY << endl;
	cout << newMinZ << " " << newMaxZ << endl;
	if ( !k_SetOutsideSelectedRegionToZero )
	{
		outfile << newNX << " " << newNY << " " << newNZ << endl;	
		outfile << newMinX << " " << newMaxX << endl;
		outfile << newMinY << " " << newMaxY << endl;
		outfile << newMinZ << " " << newMaxZ << endl;
	}
	else
	{
		outfile << nX << " " << nY << " " << nZ << endl;	
		outfile << minX << " " << maxX << endl;
		outfile << minY << " " << maxY << endl;
		outfile << minZ << " " << maxZ << endl;
	}
	
	double posX, posY, posZ;
	if (isCT)
	{
		int matidx;
		for (int iZ = 0; iZ < nZ; iZ++)
		{
			for (int iY = 0; iY < nY; iY++)
			{
				for (int iX = 0; iX < nX; iX++)
				{
					ffile >> matidx;
					posX = minX + iX * binSizeX + 0.5 * binSizeX;
					posY = minY + iY * binSizeY + 0.5 * binSizeY;
					posZ = minZ + iZ * binSizeZ + 0.5 * binSizeZ;
					if (    posX > newMinX && posX < newMaxX
					     && posY > newMinY && posY < newMaxY
					     && posZ > newMinZ && posZ < newMaxZ )
					{
						outfile << matidx << " ";
					}
					else if ( k_SetOutsideSelectedRegionToZero )
					{
						outfile << 0 << " ";
					}
				}
				if ( !k_SetOutsideSelectedRegionToZero )
				{
					if (    posY > newMinY && posY < newMaxY
				   	     && posZ > newMinZ && posZ < newMaxZ ) outfile << endl;
				}
				else
				{
					if (    posY > minY && posY < maxY
				   	     && posZ > minZ && posZ < maxZ ) outfile << endl;
				}
			}
		}
	}

	double matdensity;	// or in the case of PET: activity
	for (int iZ = 0; iZ < nZ; iZ++)
	{
		for (int iY = 0; iY < nY; iY++)
		{
			bool ok = k_SetOutsideSelectedRegionToZero;
			int iNewX(0);
			for (int iX = 0; iX < nX; iX++)
			{
				ffile >> matdensity;
				posX = minX + iX * binSizeX + 0.5 * binSizeX;
				posY = minY + iY * binSizeY + 0.5 * binSizeY;
				posZ = minZ + iZ * binSizeZ + 0.5 * binSizeZ;

				if (    posX > newMinX && posX < newMaxX
				     && posY > newMinY && posY < newMaxY
				     && posZ > newMinZ && posZ < newMaxZ )
				{
					ok = true;
					outfile << matdensity << " ";
					if ( !k_SetOutsideSelectedRegionToZero && isCT)
					{
						// if ( (iX+1)%8 == 0) outfile << endl;
						if ( (iNewX+1)%8 == 0) outfile << endl;
						iNewX++;
					}
				}
				else if ( k_SetOutsideSelectedRegionToZero )
				{
					outfile << 0.0 << " ";
					if ( isCT)
					{
						if ( (iX+1)%8 == 0) outfile << endl;
					}
				}
			}
			if (!isCT && ok)
				outfile << endl;
		}
	}
}





