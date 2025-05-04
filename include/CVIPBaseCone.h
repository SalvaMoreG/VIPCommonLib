
#ifndef __CC_VIP_BASECONE_H___
#define __CC_VIP_BASECONE_H___

#include "CVIP3Vector.h"
#include "CVIPFieldOfView.h"

#include <iostream>
#include <fstream>

class CVIPBaseCone
{
public:
							CVIPBaseCone();
 	virtual					~CVIPBaseCone();

	// copy and assign constructors
							CVIPBaseCone(const CVIPBaseCone& in_obj);
	CVIPBaseCone& 			operator= (const CVIPBaseCone& in_obj);

	// returns error
	virtual int				Set(  const C3Vector& in_position1, const double& in_e1
				   				, const C3Vector& in_position2, const double& in_e2
								, const double& in_Etot );

	inline double			GetComptonAngle() const 			{ return m_comptonAngle; }
	inline const C3Vector&	GetComptonAxisDirection() const 	{ return m_comptonAxisDirection; }
	inline const C3Vector&	GetComptonAxisOrigin() const 		{ return m_comptonAxisOrigin; }

	double					GetE1() const { return m_E1; }
	double					GetE2() const { return m_E2; }	

protected:

private:
	int 					CalculateComptonAngle(const double& in_e1, const double& in_e2, const double& in_Etot);
	int						CalculateComptonGeometrics(const C3Vector& in_position1, const C3Vector& in_position2);

	// private data
	double					m_comptonAngle;
	C3Vector				m_comptonAxisDirection;
	C3Vector				m_comptonAxisOrigin;
	
	double					m_E1;
	double					m_E2;	

	// define output stream operator
	friend std::ostream& 	operator<<(std::ostream& os, const CVIPBaseCone& in_cone);
};


#endif

