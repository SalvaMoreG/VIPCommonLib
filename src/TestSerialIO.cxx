#include "CVIPSerialIO.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// Reads in a (user given) ascii LM file
// Writes out a new binary LM file
// Reads in the same binary LM file
// Writes out a new ascii LM file

int main()
{
    CVIPSerialIO serialIO;
    
    cout << "Give input LM file" << endl;
    string fname; 
    cin >> fname;
    std::ifstream infile(fname.c_str());
    fname = "LM_binary.out";
    std::ofstream outfile(fname.c_str());

    // variables
    // unsigned long long eventID;
    unsigned long eventID; // max = 4294967295
	unsigned long long int timestamp_ps;
	std::string detString;
    double energy_MeV;
	double posX, posY, posZ;
    // float posX, posY, posZ;

    // https://www.quora.com/Does-a-text-file-have-the-same-file-size-as-a-binary-file-if-they-both-contained-the-same-integers    
    // [...] binary files are typically smaller than the same numbers in text, 
    // except in the cases when the average length of the text numbers (plus separator) is shorter 
    // than the size of the binary numbers. That typically happens if you use a larger binary format
    // than required, or if your numbers are typically small but there are a few huge outliers.
    //
    // with unsigned long long int and many doubles, the binary file will get bigger than the ascii 

    
    while ( !infile.eof() )
    {
        // serialIO.SerialRead(infile, isbinary);
        serialIO.SerialRead(infile, false, 
                            eventID, detString, timestamp_ps, energy_MeV, posX, posY, posZ); 
        if ( !infile.eof() )
        {
            serialIO.SerialWrite(outfile, true, 
                        eventID, detString, timestamp_ps, energy_MeV, posX, posY, posZ);
        }
    }
    infile.close();
    outfile.close();
    
    infile.open(fname.c_str());     // fname = "LM_binary.out";
    fname = "LM_ascii.out";
    outfile.open(fname.c_str()); 
    while ( !infile.eof() )
    {
        serialIO.SerialRead(infile, true, 
                        eventID, detString, timestamp_ps, energy_MeV, posX, posY, posZ); 
        // cout << "read binary: " << eventID << " " << detString << " " 
        //     << setprecision(15) << energy_MeV << endl;
        if ( !infile.eof() )
        {
            serialIO.SerialWrite(outfile, false, 
                        eventID, detString, timestamp_ps, energy_MeV, posX, posY, posZ);
        }
    }
    infile.close();
    outfile.close();
    
    return 0;
}


