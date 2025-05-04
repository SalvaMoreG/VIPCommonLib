#pragma once

#ifndef _VIPCOMMON_CHIT___
#define _VIPCOMMON_CHIT___

#include "CVIP3Vector.h"
#include <string>
#include <fstream>

class CHit
{
public:
						CHit();
						CHit(const C3Vector& in_position, const double& in_e);
	virtual				~CHit();

	// copy and assign constructors
						CHit(const CHit& in_obj);
	CHit& 				operator= (const CHit& in_obj);

	void 				Set(const C3Vector& in_position, const double& in_e);

	inline double				GetE() const 				{ return m_e; }
	inline const C3Vector&		GetPosition() const			{ return m_position; }

	bool				operator==(const CHit& in_obj) const; 	// ==
	bool				operator!=(const CHit& in_obj) const;	// !=

	// (De)Serialization
	void				Serialize( std::ostream& outstream ) const;
	void				Deserialize( std::istream& instream );

private:

	// define output stream operator
	friend std::ostream& 		operator<<( std::ostream& os, const CHit& in_hit);

	// data
	double 				m_e;
	C3Vector			m_position;	
};

#endif
