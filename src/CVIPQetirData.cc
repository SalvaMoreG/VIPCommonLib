
#include "CVIPQetirData.h"

#include <iostream>
#include <cstdlib>

using namespace std;

void
CVIPQetirData::Write( const C3Vector& in_pos1, const C3Vector& in_pos2
                        , const float& in_deltaT, FILE_FORMAT inFileFormat )
{
    Write( in_pos1, in_pos2, inFileFormat );
	if (inFileFormat == FFORMAT_BINARY)
	{
    	m_outputFile.write((char*) &in_deltaT, sizeof(in_deltaT));
	}
	else
	{
		m_outputFile << in_deltaT << std::endl;
	}
}

void
CVIPQetirData::Write( const C3Vector& in_pos1, const C3Vector& in_pos2
                        , FILE_FORMAT inFileFormat )
{
    if ( !m_outputFile.is_open() )
    {
        std::cout << "ERROR! Output file not open" << std::endl;
        exit(1);
    }
    float value;
    const C3Vector* position = &in_pos1;
    for (int i = 0; i < 2; i++)
    {
        if (inFileFormat == FFORMAT_BINARY)
        {
            value = position->GetX();
            m_outputFile.write((char*) &value, sizeof(value));
            value = position->GetY();
            m_outputFile.write((char*) &value, sizeof(value));
            value = position->GetZ();
            m_outputFile.write((char*) &value, sizeof(value));
        }
        else if (inFileFormat == FFORMAT_ASCII_LM)
        {
            m_outputFile << *position << " ";
            if (i == 1)
                m_outputFile << std::endl;
        }

        position = &in_pos2;
    }
}



