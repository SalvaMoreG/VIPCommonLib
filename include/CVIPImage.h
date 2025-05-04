
#ifndef __CC_VIP_VipImage_H___
#define __CC_VIP_VipImage_H___

#include <vector>
#include "CVIPFieldOfView.h"

enum FILE_FORMAT
{
	FFORMAT_NONE = 0,
	FFORMAT_ASCII_LM,
	FFORMAT_ASCII_BINNED,
	FFORMAT_BINARY
};

enum DATA_FORMAT
{
	DFORMAT_NONE = 0,
	DFORMAT_FLOAT,
	DFORMAT_UNSIGNEDINT
};

class CVipImage
{
public:
	// constructor and destructor
						CVipImage();
						CVipImage( int size );
						~CVipImage();

	// copy and assign constructors
						CVipImage(const CVipImage& in_obj);
	CVipImage& 			operator= (const CVipImage& in_obj);

	// Read/Write function
	void				Read( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
							, FILE_FORMAT in_fileformat, DATA_FORMAT in_dataformat);

	void				Write( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
							, FILE_FORMAT in_fileformat ) const;

	// Get/Set functions to access the value of a voxel
	double 				GetVoxelValue( int in_voxel ) const;
	void	 			SetVoxelValue( int in_voxel, const double& in_value );

	double				GetTotalContents() const;

	int					GetSize() const { return m_size; }

	// operators
	void				Multiply( const double& in_factor );
	CVipImage 			operator+(const CVipImage&) const;       // operator+()

protected:

private:
    // private functions
    void                ReadBinaryImageFile(const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
											, DATA_FORMAT in_format);
    void                ReadAsciiImageFile(const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView
											, bool is_ListMode );

	void				OutputImageBinary( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView ) const;
	void				OutputImageBinnedAscii( const std::string& in_filename, const CVIPFieldOfView& in_fieldOfView ) const;

protected:
    // data
    
	int					m_size;
	std::vector<double> m_imgdata;

};

#endif

