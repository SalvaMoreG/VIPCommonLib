    
#ifndef __CC_VIP_BASELOR_H___
#define __CC_VIP_BASELOR_H___

#include "CVIP3Vector.h"
#include "CVIPFieldOfView.h"

#include <vector>

class CVIPBaseLor
{
public: 
							CVIPBaseLor();
							CVIPBaseLor(  const C3Vector& in_position1, const double& in_e1
				   			    , const C3Vector& in_position2, const double& in_e2 );
	virtual					~CVIPBaseLor();

	virtual int				Set(  const C3Vector& in_position1, const double& in_e1
				   			    , const C3Vector& in_position2, const double& in_e2 );

	inline const C3Vector& 	GetPosition1() const { return m_position1; }
	inline const C3Vector& 	GetPosition2() const { return m_position2; }
	
	double					GetE1() const { return m_E1; }
	double					GetE2() const { return m_E2; }	

	void					FindFovBins( const CVIPFieldOfView& in_fieldOfView
								, std::vector<unsigned int>& io_fovidxVec
								, std::vector<float>* io_fovweightsVec = 0 ) const;
//								, std::vector<unsigned int>* io_fovweightsVec = 0 ) const;

	enum LOROrientation
	{
		LORDIR_NONE = 0,
		LORDIR_X, 
		LORDIR_Y, 
		LORDIR_Z
	};
	inline LOROrientation	GetOrientation() const { return m_lorOrientation; }

private:
	// aux
	double 					GetMinFactorPromille( const CVIPFieldOfView& in_fieldOfView ) const;

	// define output stream operator
	friend std::ostream& 		operator<<( std::ostream& os, const CVIPBaseLor& in_lor);
	
	// private data
	C3Vector				m_position1;
	C3Vector				m_position2;
	
	double					m_E1;
	double					m_E2;

	LOROrientation			m_lorOrientation;
};

// ==========================================================================================


#endif

