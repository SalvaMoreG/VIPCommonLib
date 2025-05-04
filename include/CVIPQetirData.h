
#pragma once

#ifndef __CVIP_QetirData_H___
#define __CVIP_QetirData_H___

#include "CVIPHit.h"
#include "CVIP3Matrix.h"
#include "CVIPImage.h"

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

class CVIPQetirData
{
public:
					CVIPQetirData()
					{
					}

	virtual			~CVIPQetirData()
					{
						Finalize();
					}

    void            OpenOutputFile( const std::string& in_filename )
                    {
                        Initialize( in_filename );
                    }
    void            CloseOutputFile()
                    {
                        Finalize();
                    }

	void			Write( const C3Vector& in_pos1, const C3Vector& in_pos2
                        , const float& in_deltaT, FILE_FORMAT inFileFormat );

	void			Write( const C3Vector& in_pos1, const C3Vector& in_pos2
                        , FILE_FORMAT inFileFormat = FFORMAT_BINARY );

private:
	//	Initialize
	void			Initialize( const std::string& in_filename )
					{
						m_outputFile.open( in_filename.c_str(), std::ios::binary | std::ios::out );
						assert( m_outputFile.is_open() );
					}

	void			Finalize()
					{
						m_outputFile.close();
					}

	// data
	std::ofstream	m_outputFile;
	C3Matrix 		m_matrix;
};

#endif


