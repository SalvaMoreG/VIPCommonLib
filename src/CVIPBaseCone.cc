
#include "CVIPBaseCone.h"

using namespace std;

ostream& operator<<(ostream& os, const CVIPBaseCone& in_cone)
{
	os << "Angle: " << in_cone.GetComptonAngle() << " = " << in_cone.GetComptonAngle()*180/kPI << " degrees"
       << ", Axis: " << in_cone.GetComptonAxisDirection()
       << ", Apex: " << in_cone.GetComptonAxisOrigin() 
	   << ", E1: " << in_cone.GetE1() 
	   << endl;

    return os;
}


CVIPBaseCone::CVIPBaseCone()
	: m_comptonAngle(0.0)
	, m_comptonAxisDirection(0.0, 0.0, 0.0)
	, m_comptonAxisOrigin(0.0, 0.0, 0.0)
	, m_E1(0.0)
	, m_E2(0.0)
{
}

CVIPBaseCone::~CVIPBaseCone()
{
}

// copy constructor
CVIPBaseCone::CVIPBaseCone(const CVIPBaseCone& in_obj)
{
	*this = in_obj;
}

// assignment operator
CVIPBaseCone&
CVIPBaseCone::operator= (const CVIPBaseCone& in_obj)
{
	if (this != &in_obj)
	{
		m_comptonAngle = in_obj.GetComptonAngle();
		m_comptonAxisDirection = in_obj.GetComptonAxisDirection();
		m_comptonAxisOrigin = in_obj.GetComptonAxisOrigin();
		m_E1 = in_obj.GetE1();
		m_E2 = in_obj.GetE2();
	}
	return *this;
}

int
CVIPBaseCone::Set(  const C3Vector& in_position1, const double& in_e1
        , const C3Vector& in_position2, const double& in_e2
        , const double& in_Etot )
{
    int error = CalculateComptonAngle(in_e1, in_e2, in_Etot);
	if (!error)
	{
		CalculateComptonGeometrics(in_position1, in_position2);
	
		m_E1 = in_e1;
		m_E2 = in_e2;
	}
	
	return error;
}

int
CVIPBaseCone::CalculateComptonAngle(const double& in_e1, const double& in_e2, const double& in_Etot)
{
    // Etot
	// double Etot = CUserParameters::Instance()->GetGammaSourceEnergy();
	//
	double cosTh;
	if (in_Etot <= 0.0)
		cosTh = 1.0 - mass_electron_keV*((1/(in_e2)) - (1/(in_e1+in_e2)));
	else 
		cosTh = 1.0 - mass_electron_keV*((1/(in_Etot - in_e1)) - (1/(in_Etot)));

	if (fabs(cosTh) > 1.0)
	{
		return 1;	// ERROR
	}

	m_comptonAngle = acos(cosTh);

	return 0;
}

int
CVIPBaseCone::CalculateComptonGeometrics(const C3Vector& in_position1, const C3Vector& in_position2)
{
    m_comptonAxisDirection = in_position1 - in_position2;
	//	cout << "pos1: " << in_position1 << " pos2: " << in_position2 << " so AXIS: " << m_comptonAxisDirection << endl;
	double len = m_comptonAxisDirection.GetLength();
	if (len > 0.0)
	{
		m_comptonAxisDirection = m_comptonAxisDirection * ( 1.0 / len );
	}

	m_comptonAxisOrigin = in_position1;
	return 0;
}




