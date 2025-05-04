#include "CVIPSerialIO.h"

#include "CVIPUtils.h"

#include <iostream>
#include <cassert>
#include <string>
#include <iomanip>

using namespace std;

template <typename T, typename... Types>
void CVIPSerialIO::SerialRead(std::ifstream& infile, bool is_binary
        , T& var1, Types&... var2)
{
    if (is_binary)
        VIPUtils::DoReadBinary(infile, var1);
    else 
        infile >> var1;
    //  cout << "read var1: " << var1 << " ";
	if ( !infile.eof() )
	{
        // call SerialRead with next variable in list
		SerialRead(infile, is_binary, var2...);
	}    
}

template <typename T, typename... Types>
void CVIPSerialIO::SerialWrite(std::ofstream& outfile, bool is_binary
        , const T& var1, const Types&... var2)
{
    if (is_binary)
        VIPUtils::DoWriteBinary(outfile, var1);
    else
        outfile << setprecision(15) << var1 << " ";

	// call SerialWrite with next variable in list
	SerialWrite(outfile, is_binary, var2...);
}


