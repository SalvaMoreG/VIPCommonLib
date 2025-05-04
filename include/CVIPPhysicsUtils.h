
#ifndef _C_PHYSICS_UTILS__H__
#define _C_PHYSICS_UTILS__H__

#include "CVIP3Vector.h"
#include "VIPconstants.h"
#include <vector>

// TODO:
// Make this is a Namespace (but with data???) or a real class....

class CPhysicsUtils
{
public:
							CPhysicsUtils();
							~CPhysicsUtils();

	static double 			CalculateKleinNishinaProb( const double& in_Etot_keV, const double& in_costh );

    double					GetP3ScattererComptonProbability( const double& in_E_keV, const double& pathlength ) const;
	double					GetP6AbsorberPhotoElectricProbability( const double& in_E_keV, const double& pathlength ) const;
    double					GetP5ScattererEscapeProbability( const double& in_E_keV, const double& pathlength ) const;

	double					GetTotalAttenuation( const double& in_energy_keV, bool isCdTe ) const;

private: 

	void					LookUpTable( bool isCdTe );
	void					InterpolateEnergies( bool isCdTe, const double& in_E_keV
											  , int& io_index, double& io_prevE, double& io_nextE ) const;

	std::vector<double>		m_CdTeEnergies;
	std::vector<double>		m_CdTePhotoElAtt;
	std::vector<double>		m_CdTeTotWiAtt;

	std::vector<double>		m_SiEnergies;
	std::vector<double>		m_SiComptonAtt;
	std::vector<double>		m_SiTotWiAtt;

	double 					m_maxP6;
	double 					m_maxP3;
};

// global functions

// ------------------------------------------------------------------------------------------

inline double ComptonCosThFromE1(const double& in_E1, const double& in_Etot)
{
	if (in_E1 == 0) return 0;

	double term = in_Etot - in_E1;
	double cosTh = 1 - mass_electron_keV * ( 1./term - 1/in_Etot);
	return cosTh;
}

// ------------------------------------------------------------------------------------------

inline double ComptonE1FromCosTh(const double& in_cosTh, const double& in_Etot)
{
	double term = in_Etot * (1 - in_cosTh) + mass_electron_keV;
	double E1 = in_Etot * (1.0 - mass_electron_keV / term);
	return E1;
}


#endif

