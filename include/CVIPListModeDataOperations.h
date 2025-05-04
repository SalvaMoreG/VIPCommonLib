
#pragma once

#ifndef __CVIP_ListModeDataOperations_H___
#define __CVIP_ListModeDataOperations_H___

#include "CVIPHit.h"
#include "CVIP3Matrix.h"
#include "CVIPImage.h"

#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

class CVIPListModeDataOperations
{
public:
					CVIPListModeDataOperations() {}

	virtual			~CVIPListModeDataOperations() {}

	void			TangentialBlurring( C3Vector& io_position );
	void			TangentialBlurringPEM( C3Vector& io_position ) const; 

private:
	// data
	C3Matrix 		m_matrix;	// only created once, initialized each time TB function is called.
};

#endif


