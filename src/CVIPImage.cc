
#include "CVIPImage.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

CVipImage::CVipImage()
    : m_size(0)
	, m_imgdata(0)    
{
}

CVipImage::CVipImage( int in_size )
    : m_size(in_size)
	, m_imgdata(0)
{
	m_imgdata.resize(m_size, 0.0);
}

// ****************************************************************

CVipImage::~CVipImage()
{
}

// ****************************************************************

// copy constructor
CVipImage::CVipImage(const CVipImage& in_obj)
    : m_size(0)
	, m_imgdata(0)
{
	*this = in_obj;
}

// ****************************************************************

// assign operator
CVipImage&
CVipImage::operator= (const CVipImage& in_obj)
{
	if (this != &in_obj)
	{
		m_size = in_obj.GetSize();
		
		m_imgdata.resize(m_size, 0);
	
		for (int ivoxel = 0; ivoxel < m_size; ivoxel++)
		{
			m_imgdata[ivoxel] = in_obj.GetVoxelValue( ivoxel );
		}
	}
	return *this;
}

// ****************************************************************

void
CVipImage::Read( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
					, FILE_FORMAT in_fileformat, DATA_FORMAT in_dataformat)
{
	// m_imgdata is the array that contains all the values for all voxels in the image
	m_size = in_fieldOfView.GetNumberOfVoxels();
	m_imgdata.resize(m_size, 0);

	// Read file
	if (in_fileformat == FFORMAT_ASCII_LM)
	{
		ReadAsciiImageFile( in_filename, in_fieldOfView, true );
	}
	else if (in_fileformat == FFORMAT_ASCII_BINNED)
	{
		ReadAsciiImageFile( in_filename, in_fieldOfView, false );
	}
	else if (in_fileformat == FFORMAT_BINARY)
	{
		ReadBinaryImageFile( in_filename, in_fieldOfView, in_dataformat );
	}
	else
	{
		cout << "ERROR. Wrong file format." << endl;
		exit(2);
	}
}

// ****************************************************************

void
CVipImage::Write( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
				, FILE_FORMAT in_fileformat ) const
{
	if ( in_fileformat == FFORMAT_BINARY)
	{
		OutputImageBinary( in_filename, in_fieldOfView );
	}
	else if ( in_fileformat != FFORMAT_NONE )
	{
		OutputImageBinnedAscii( in_filename, in_fieldOfView );
	}
}


// ****************************************************************

double
CVipImage::GetVoxelValue( int in_voxel ) const
{
	if (in_voxel < 0 || in_voxel >= m_size)
	{
		cout << "in_voxel: " << in_voxel << " whereas size: " << m_size << endl;
	}
	assert (in_voxel >= 0 && in_voxel < m_size);
	return m_imgdata[in_voxel];
}

// ****************************************************************

void
CVipImage::SetVoxelValue( int in_voxel, const double& in_value )
{
	if (in_voxel < 0 || in_voxel >= m_size)
	{
		cout << "CIm:SetVoxelValue, in_voxel: " << in_voxel << " size: " << m_size << endl;
	}
	assert (in_voxel >= 0 && in_voxel < m_size);
	m_imgdata[in_voxel] = in_value;
}

// ****************************************************************

double
CVipImage::GetTotalContents() const
{
	double totContents = 0.0;
	for (int ivoxel = 0; ivoxel < m_size; ivoxel++)
	{
		totContents += m_imgdata[ivoxel];
	}
	return totContents;
}

// ****************************************************************

void
CVipImage::Multiply( const double& in_factor )
{
	for (int ivoxel = 0; ivoxel < m_size; ivoxel++)
	{
		m_imgdata[ivoxel] = m_imgdata[ivoxel] * in_factor;
	}
}

// ****************************************************************

// + (plus) operator
CVipImage
CVipImage::operator+ (const CVipImage& in_obj) const
{
	assert ( m_size == in_obj.GetSize() );

    CVipImage newImage( m_size );
    double new_value = 0.0;
	for (int ivoxel = 0; ivoxel < m_size; ivoxel++)
	{
		new_value = m_imgdata[ivoxel] + in_obj.GetVoxelValue( ivoxel );
		newImage.SetVoxelValue(ivoxel, new_value);
	}

    return newImage;
}

// ****************************************************************

void
CVipImage::ReadBinaryImageFile(const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView, DATA_FORMAT in_format)
{
	ifstream binaryfile( in_filename.c_str(), ios::binary | ios::in );
	if (!binaryfile.is_open())
	{
		cout << "ERROR. Cannot open file: " << in_filename << endl;
		exit(1);
	}

	int tmpint;     // TODO: shouldn't this be a "unsigned int"??????
	float tmpfloat;
	int ivoxel = 0;

    int nbinsX = in_fieldOfView.GetNumVoxelsX();
    int nbinsY = in_fieldOfView.GetNumVoxelsY();
    int nbinsZ = in_fieldOfView.GetNumVoxelsZ();

    // loop over all bins in the AMIDE file
	for (int kz=0; kz < nbinsZ; kz++)
	{
		for (int jy=0; jy < nbinsY; jy++)
		{
			for (int ix=0; ix < nbinsX; ix++)
			{
				ivoxel = in_fieldOfView.GetVoxelIndex( ix, jy, kz );

				if (in_format == DFORMAT_UNSIGNEDINT)
				{
					// Read voxel value from file
					binaryfile.read((char*) &tmpint, sizeof(tmpint));

					// Fill array
					if ( ivoxel >= 0 && ivoxel < m_size && tmpint > 0.0 )
					{
						m_imgdata[ivoxel] = (float) tmpint;
					}
				}
				else if (in_format == DFORMAT_FLOAT)
				{
					// Read voxel value from file
					binaryfile.read((char*) &tmpfloat, sizeof(tmpfloat));

					// Fill array
					if ( ivoxel >= 0 && ivoxel < m_size && tmpfloat > 0.0 )
					{
						m_imgdata[ivoxel] = tmpfloat;
					}
				}
			}
		}
	}
	binaryfile.close();
}

// ****************************************************************

void
CVipImage::ReadAsciiImageFile(const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView, bool is_ListMode )
{
	std::ifstream datfile (in_filename.c_str(), ios::in);
	if (!datfile.is_open())
	{
		cout << "ERROR. ASCII FILE NOT FOUND:" << in_filename << endl;
		exit(1);
	}

    int nbinsZ = in_fieldOfView.GetNumVoxelsZ();

	int idum, idum2, kz;
	double x, y, z;
    int ivoxel = 0;
	int nOutside = 0;
	double weight = 1.0;

	while (!datfile.eof())
	{
		if (is_ListMode)
			datfile >> x >> y >> z;
		else 
			datfile >> x >> y >> z >> weight;

		if (!datfile.eof())
		{
			if ( in_fieldOfView.IsInBounds( C3Vector(x, y, z) ) == false )	// if position is NOT inside bounds of FOV
			{
				nOutside++;
			}
			else
			{
				ivoxel = in_fieldOfView.GetVoxelIndex( C3Vector(x, y, z) );
				if (ivoxel >= 0 && ivoxel < m_size)
				{
					m_imgdata[ivoxel] += weight;   
						// Listmode: add an other count (=event) to this voxel
						// Binmode: add bin contents to this voxel
				}
			}
		}
	}
    datfile.close();
}

// ****************************************************************

void
CVipImage::OutputImageBinary( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView ) const
{
	ofstream binaryfile(in_filename.c_str(), ios::binary | ios::out );
	if (!binaryfile.is_open())
	{
		cout << "ERROR. Cannot open file: " << in_filename << endl;
		exit(1);
	}

    int nbinsX = in_fieldOfView.GetNumVoxelsX();
    int nbinsY = in_fieldOfView.GetNumVoxelsY();
    int nbinsZ = in_fieldOfView.GetNumVoxelsZ();

	int ivoxel;
	for (int kz=0; kz < nbinsZ; kz++)
	{
		for (int jy=0; jy < nbinsY; jy++)
		{
			for (int ix=0; ix < nbinsX; ix++)
			{
				ivoxel = in_fieldOfView.GetVoxelIndex( ix, jy, kz );

				float tmp = m_imgdata[ivoxel];
				binaryfile.write((char*) &tmp, sizeof(tmp));
			}
		}
	}
	binaryfile.close();
}

// ****************************************************************

void
CVipImage::OutputImageBinnedAscii( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView ) const
{
	ofstream ofile(in_filename.c_str() );
	if (!ofile.is_open())
	{
		cout << "ERROR. Cannot open file: " << in_filename << endl;
		exit(1);
	}

    int nbinsX = in_fieldOfView.GetNumVoxelsX();
    int nbinsY = in_fieldOfView.GetNumVoxelsY();
    int nbinsZ = in_fieldOfView.GetNumVoxelsZ();

	float weight;
	int ivoxel;
	double x, y, z;
	for (int kz=0; kz < nbinsZ; kz++)
	{
		for (int jy=0; jy < nbinsY; jy++)
		{
			for (int ix=0; ix < nbinsX; ix++)
			{
				ivoxel = in_fieldOfView.GetVoxelIndex( ix, jy, kz );
				weight = m_imgdata[ivoxel];
				{
					in_fieldOfView.GetVoxelCentre( ix, jy, kz, x, y, z );
					ofile << x << " " << y << " " << z << "  " << weight << endl;
				}
			}
		}
	}
	ofile.close();
}

// *******************************************************************************






