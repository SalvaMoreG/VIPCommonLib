
#include "CVIPImageKernel.h"

#include <iostream>
#include <cmath>

#include "VIPconstants.h"
#include "CVIP3Vector.h"
#include "CVIP3Matrix.h"
#include "CVIPUtils.h"


/*
// From: "http://en.wikipedia.org/wiki/Kernel_(image_processing)"
//
	for each image row in output image:
		for each pixel in image row:
	
			set accumulator to zero
	
			for each kernel row in kernel:
				for each element in kernel row:
	
					if element position  corresponding* to pixel position then
						multiply element value  corresponding* to pixel value
						add result to accumulator
					endif
	
			set output image pixel to accumulator
*/

using namespace std;

// ======================================================================================================

CVIPGaussKernel::CVIPGaussKernel( int in_size, int in_dimension, const double& in_sigma )
	: m_data(NULL)
	, m_size(in_size)
	, m_dimension(in_dimension)
	, m_sigma(in_sigma)		
	, m_datasize(0)
	, m_amplitude(0)
{
	// "normalize" Gauss sigma on "standard matrix" of 1x1x1 mm^3 voxels
	// So, for instance, with a kernel of size 3, we could use 2*sigma (= 95.4%) = 1.5 ---> sigma = 0.75

	Initialize();
}

// ======================================================================================================

CVIPGaussKernel::~CVIPGaussKernel()
{
	delete [] m_data;
}

// ======================================================================================================

double
CVIPGaussKernel::Convolute( const CVIPVector& in_image, const CVIPFieldOfView& in_fov, int in_index )
{
	return DoConvolute( in_image, in_fov, in_index );
}

// ======================================================================================================

void
CVIPGaussKernel::PrintGauss() const
{
	int izmax = (m_dimension == 2) ? 0 : m_size-1;
	int index_kernel;
	for (int istepZ = 0; istepZ <= izmax; istepZ++)
	{
		for (int istepY = 0; istepY <= m_size-1; istepY++)
		{
			for (int istepX = 0; istepX <= m_size-1; istepX++)
			{
				index_kernel = GetIndex(istepX, istepY, istepZ);
				cout << m_data[index_kernel] << "     ";
			}
			cout << endl;
		}
	}
	if (m_amplitude > 0)
	{
		cout << "SUM: " << 1./m_amplitude << " with amplitude: " << m_amplitude << endl;
	}
	else
		cout << "WRONG amplitude " << m_amplitude << endl;
}

// ======================================================================================================

void
CVIPGaussKernel::Initialize()
{
	m_datasize = m_size;
	for (int i = 1; i < m_dimension; i++)
		m_datasize = m_datasize * m_size;

	m_data = new double [m_datasize];
	
	// Gaussian parameters
	C3Vector mean(0.0, 0.0, 0.0);

	// fill data with Gaussian values...
	double g = VIPUtils::Gauss( mean, m_sigma, C3Vector(0.0, 0.0, 0.0), 2);
	double A = 4./g;	// whatever....

	// cout << " sigma: " << m_sigma << endl;
	// cout << " A: " << A << " g: " << g << endl;

	m_amplitude = 0.0;
	int index_kernel;
	int offset = -1 * (m_size>>1);
	int izmax = (m_dimension == 2) ? 0 : m_size-1;

	double x0 = mean.GetX() + offset;
	double y0 = mean.GetY() + offset;
	double z0 = mean.GetZ();
	if (m_dimension > 2) z0 += offset;
	double x, y, z;

	for (int istepZ = 0; istepZ <= izmax; istepZ++)
	{
		z = z0 + istepZ;
			// Gauss is made on a matrix with pixel-size 1x1x1 (independent of voxelsize of image-FOV)

		for (int istepY = 0; istepY <= m_size-1; istepY++)
		{
			y = y0 + istepY;
				// Gauss is made on a matrix with pixel-size 1x1x1 (independent of voxelsize of image-FOV)
			
			for (int istepX = 0; istepX <= m_size-1; istepX++)
			{
				x = x0 + istepX;
					// Gauss is made on a matrix with pixel-size 1x1x1 (independent of voxelsize of image-FOV)

				g = A * VIPUtils::Gauss( mean, m_sigma, C3Vector(x, y, z), m_dimension);

				/*
				if (y == 0 && z == 0)
					cout << " position: " << C3Vector(x, y, z) << " AND g:" << g << endl;
				*/

				index_kernel = GetIndex(istepX, istepY, istepZ);
				m_data[index_kernel] = g;	

				m_amplitude += g;
			}
			// cout << endl;
		}
		// cout << "		---         " << endl;
	}
	if (m_amplitude != 0.0)
		m_amplitude = 1./m_amplitude;
}

// ======================================================================================================

double
CVIPGaussKernel::DoConvolute( const CVIPVector& in_image, const CVIPFieldOfView& in_fov, int in_index )
{
	int ix0, iy0, iz0;
	in_fov.GetVoxelIndices( in_index, ix0, iy0, iz0);

	int index_img;
	int index_kernel;
	int offset = -1 * (m_size>>1);
	int izmax = (m_dimension == 2) ? 0 : m_size-1;

	ix0 += offset;
	iy0 += offset;
	if (m_dimension > 2) iz0 += offset;

	int ix, iy, iz;
	double sum = 0.0;

bool dbg = false;
if (dbg) 
{
cout << endl;
cout << "in_index: " << in_index << " so center xyz: " << ix0 << " " << iy0 << " " << iz0 << endl;
cout << "izmax: " << izmax << " m_size: " << m_size << " offset: " << offset << endl;
}	

	for (int istepZ = 0; istepZ <= izmax; istepZ++)
	{
		iz = iz0 + istepZ;
if (dbg) cout << "iz: " << iz << endl;
		if (iz >= 0 && iz < in_fov.GetNumVoxelsZ())
		{
			for (int istepY = 0; istepY <= m_size-1; istepY++)
			{
				iy = iy0 + istepY;
				if (iy >= 0 && iy < in_fov.GetNumVoxelsY())
				{
					for (int istepX = 0; istepX <= m_size-1; istepX++)
					{
						ix = ix0 + istepX;
						if (ix >= 0 && ix < in_fov.GetNumVoxelsX())
						{
							index_img = in_fov.GetVoxelIndex(ix, iy, iz);
							index_kernel = GetIndex(istepX, istepY, istepZ);
	
							sum += in_image[index_img]*m_data[index_kernel]*m_amplitude;
if (dbg) cout << "ixy: " << ix << " " << iy << " index_kernel: " << index_kernel << " index_img: " << index_img 
	 << " img: " << in_image[index_img] << " kernel*A: " << m_data[index_kernel]*m_amplitude << " sum: " << sum << endl;
						}
					}
				}
			}
		}
	}
	return sum;
}

// ======================================================================================================

int
CVIPGaussKernel::GetIndex(int in_ix, int in_iy, int in_iz) const
{
	return in_iz*m_size*m_size + in_iy*m_size + in_ix;
}

// ======================================================================================================

