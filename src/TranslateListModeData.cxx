
// In this case, we assume there is no intersection, i.e. no "delta t"...

#include "CVIPQetirData.h"
#include "CVIPListModeDataOperations.h"

#include <string>
#include <iostream>
#include <fstream>

#include <cassert>
#include <cstdlib>

using namespace std;

bool CheckNumberOfColumns( ifstream& in_file, int in_nrReqColumns, const string& in_filename );

int
main(int argc, char* argv [])
{
	cout << "give 'PET_img.out' filename" << endl;
	string infilename;
	cin >> infilename;
	cout << "LIST MODE INPUT FILE: " << infilename << endl;

	ifstream infile( infilename.c_str() );
	if ( !infile.is_open() )
	{
		cout << "file not open: " << infilename << endl;
		exit(1);
	}

	bool doQETIR( true );
    cout << "Output file is QETIR data (1/0)?" << endl;
    cin >> doQETIR;
	cout << "OUTPUT FILE IS QETIR DATA: " << doQETIR << endl;

	bool doAscii( false );
	if (doQETIR) 
	{
		cout << "Make QETIR input file in ASCII format? (testing only) 1/0" << endl;
		cin >> doAscii;
		cout << "MAKE QETIR INPUT FILE IN ASCII FORMAT: " << doAscii << endl;
	}

	bool doTimeOfFlight( false );
	cout << "Also translate TOF? (LMinput file must have delta time field) 1/0" << endl;
	cin >> doTimeOfFlight;
	cout << "INCLUDE TOF (LM in/out put have TOF field): " << doTimeOfFlight << endl;

    bool doTangentialBlurring( true );
    cout << "Apply tangential blurring (instead of using voxel centers?) (1/0)" << endl;
    cin >> doTangentialBlurring;
	cout << "DO TANGENTIAL BLURRING: " << doTangentialBlurring << endl;

	bool hasAddedSrcPos( false );
	cout << "Does input file have 4 additional fields containing evtNr plus real source position? (1/0)" << endl;
	cin >> hasAddedSrcPos;
	cout << "INPUT FILE HAS ADDITIONAL INFO ABOUT SRC POSITION: " << hasAddedSrcPos << endl;

	// bool hasAddedRandomFlag( false );
	// cout << "Does input file have 1 additional field containing random-flag? (1/0)" << endl;
	// cin >> hasAddedRandomFlag;
	// cout << "INPUT FILE HAS ADDITIONAL RANDOM-FLAG: " << hasAddedRandomFlag << endl;

	double deltaZ(0.0);
	cout << "Shift z coordinate (delta-z / 0)" << endl;
	cin >> deltaZ;
	cout << "z-coordinates shifted with distance: " << deltaZ << endl;

	double deltaY(0.0);
	cout << "Shift y coordinate (delta-y / 0)" << endl;
	cin >> deltaY;
	cout << "y-coordinates shifted with distance: " << deltaY << endl;

	double deltaX(0.0);
	cout << "Shift x coordinate (delta-x / 0)" << endl;
	cin >> deltaX;
	cout << "x-coordinates shifted with distance: " << deltaX << endl;

	// if (hasAddedSrcPos && hasAddedRandomFlag)
	// {
	// 	cout << "ERROR! THIS IS TOO COMPLICATED, I QUIT!" << endl;
	// 	cout << "INPUT FILE HAS ADDITIONAL INFO ABOUT SRC POSITION: " << hasAddedSrcPos << endl;
	// 	cout << "INPUT FILE HAS ADDITIONAL RANDOM-FLAG: " << hasAddedRandomFlag << endl;
	// 	return -1;
	// }

	bool isPEM( false );
	if (doTangentialBlurring)
	{
		cout << "Is this PEM data? (1/0)" << endl;
		cin >> isPEM;
		if (isPEM)
			cout << "THIS IS PEM DATA" << endl;
		else 
			cout << "THIS IS PET DATA" << endl;
	}


	if (hasAddedSrcPos && doQETIR) 
	{
		cout << "ERROR! There is no point in having addSrcPos, when it is ignored anyway for QETIR input data" << endl;
		return -1;
	}

	if (hasAddedSrcPos && !doTangentialBlurring) 
	{
		cout << "WARNING! Output file will be identical to input file (VIP data file with NO tangential blurring" << endl;
	}

	string outfilename( "out_" );
	outfilename += infilename;
    if (doAscii)
        outfilename += "_ascii";

    if (doTangentialBlurring)
        outfilename += "_tangblur";

    if (doQETIR)
	{
		if (doTimeOfFlight)
        	outfilename += ".lmdT_NEW";
		else
			outfilename += ".lmd_NEW";
	}
    else
	{
		if (doTimeOfFlight)
        	outfilename += ".vipT_NEW";
		else
			outfilename += ".vip_NEW";
	}

	CVIPQetirData qetirData;
	ofstream outputFile;
    if (doQETIR)
    {
        qetirData.OpenOutputFile( outfilename );
    }
    else
    {
		outputFile.open( outfilename.c_str() );
        assert( outputFile.is_open() );
    }

	int iline = 0, ievent(0);
	C3Vector pos1, pos2;
	double x, y, z, E1, E2, dT;
	CVIPListModeDataOperations lmDataOperations;
	// int randomflag;
	while ( !infile.eof() )
	{
		// if (iline == 0)
		if (false)
		{
			int nrReqColumns = 8;
			if (doTimeOfFlight)
				nrReqColumns++;
			if (hasAddedSrcPos)
				nrReqColumns += 4;
			// if (hasAddedRandomFlag)
			// 	nrReqColumns += 1;
			CheckNumberOfColumns( infile, nrReqColumns, infilename );
		}

		infile >> z >> y >> x >> E1;
		if (deltaZ != 0)
			z += deltaZ;
		if (deltaY != 0)
			y += deltaY;
		if (deltaX != 0)
			x += deltaX;
        pos1.Set(x, y, z);

		infile >> z >> y >> x >> E2;
		if (deltaZ != 0)
			z += deltaZ;
		if (deltaY != 0)
			y += deltaY;
		if (deltaX != 0)
			x += deltaX;
        pos2.Set(x, y, z);
	
		if (doTimeOfFlight)
		{
			infile >> dT;
		}

		// if (hasAddedSrcPos)
		// {
		// 	infile >> ievent >> x >> y >> z;
		// }

		// if (hasAddedRandomFlag)
		// {
		// 	infile >> randomflag;
		// }

		infile.ignore(1024, '\n');

		if (doTangentialBlurring)
		{
			if (isPEM) 
			{
            	lmDataOperations.TangentialBlurringPEM( pos1 );
            	lmDataOperations.TangentialBlurringPEM( pos2 );
			}
			else
			{
            	lmDataOperations.TangentialBlurring( pos1 );
            	lmDataOperations.TangentialBlurring( pos2 );
			}
		}

		if ( !infile.eof() )
		{
            if (doQETIR)
            {
				FILE_FORMAT fileFormat = (doAscii) ? FFORMAT_ASCII_LM : FFORMAT_BINARY;
				if (doTimeOfFlight)
					qetirData.Write( pos1, pos2, dT, fileFormat );
				else
                	qetirData.Write( pos1, pos2, fileFormat );
            }
            else
            {
                outputFile << pos1.GetZ() << " " << pos1.GetY() << " " << pos1.GetX() << " " << E1 << " "
                           << pos2.GetZ() << " " << pos2.GetY() << " " << pos2.GetX() << " " << E2 << " ";
				if (doTimeOfFlight)
				{
					outputFile << dT << " ";
				}
				if (hasAddedSrcPos)
				{
					outputFile << ievent << " " << x << " " << y << " " << z;
				}
				outputFile << endl;
            }

			if (iline%1000000 == 0)
			{
				cout << "[" << iline << "]: " << pos1 << " " << pos2 << " ";
				if (doTimeOfFlight)
				{
					cout << dT << " ";
				}
				cout << endl;
			}
			iline++;
		}
	}

    if (doQETIR)
    {
        qetirData.CloseOutputFile();
    }
    else
    {
        outputFile.close();
    }
}


// ==============================================

// TODO: this function is very useful! Put it in VIPCOMMONLIB

bool CheckNumberOfColumns( ifstream& in_file, int in_nrReqColumns, const string& in_filename )
{
	int curpos = in_file.tellg();
	char dumline[256];
	in_file.getline (dumline,256);
	string strline( dumline );
	// cout << "strline: " << strline << endl;
	std::istringstream iss(strline);
	int columns = 0;
	do
	{
		std::string sub;
		iss >> sub;
		if (sub.length())
			++columns;
	}
	while(iss);
	in_file.seekg(curpos);
	cout << "FILE: " << in_filename << " #columns: " << columns << endl;

	if ( columns != in_nrReqColumns )
	{
		cout << "ERROR, wrong number of columns in FILE: " << in_filename << endl;
		cout << "first line of file: " << endl;
		cout << strline << endl;
		cout << "required #columns: " << in_nrReqColumns << endl;
		cout << "#columns in file: " << columns << endl;
	}
	assert( columns == in_nrReqColumns );

	return ( columns == in_nrReqColumns );
}











