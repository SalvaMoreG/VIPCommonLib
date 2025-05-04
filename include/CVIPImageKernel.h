
#pragma once

#ifndef _VIPCOMMON_IMAGEKERNEL_H__
#define _VIPCOMMON_IMAGEKERNEL_H__

#include "VIPconstants.h"
#include "CVIPVector.h"
#include "CVIPFieldOfView.h"

class CVIPImageKernel
{
public:
	// abstract class, never to be instantiated and no member data, so needs no constructor 
	/* 
					CVIPImageKernel();
	virtual			~CVIPImageKernel();
	*/
	double virtual	Convolute( const CVIPVector& in_image, const CVIPFieldOfView& in_fov, int image_index ) = 0;
};

// ======================================================================================================

class CVIPGaussKernel : public CVIPImageKernel
{
public: 
					CVIPGaussKernel( int size, int dimension, const double& sigma );
	virtual			~CVIPGaussKernel();

	double			Convolute( const CVIPVector& in_image, const CVIPFieldOfView& in_fov, int image_index );
	void			PrintGauss() const; 	// for debugging

private:
	
	void			Initialize();

	double			DoConvolute( const CVIPVector& in_image, const CVIPFieldOfView& in_fov, int in_index );
	int				GetIndex(int in_ix, int in_iy, int in_iz) const;

	// data
	double*			m_data;
	int 			m_size;
	int 			m_dimension;
	double			m_sigma;
	int 			m_datasize;
	double			m_amplitude;
};

// ======================================================================================================

	// double  GaussianBlur( const CVIPFieldOfView& , CVector ----> Not in VIP yet....--> CVIPVector, 
	// 		int current_index )
	// {
		/*
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
		/*
		in_fov.GetVoxelIndices( current_index, curr_x, curr_y, curr_z);

		for (int xstep = -1; xstep <= 1; xstep++)
		{
			ix += xstep;
			if (ix >= 0)
			{
				for (int ystep = -1; ystep <= 1; ystep++)
				{
					iy += ystep;
					if (iy >= 0)
					{
					}
				}
			}
		}
		*/
	// }

#endif
