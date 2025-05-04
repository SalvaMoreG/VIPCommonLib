
#pragma once

#ifndef __CVIPSEMICONDUCTORENERGY_H__
#define __CVIPSEMICONDUCTORENERGY_H__

// ****************************************
// Smearing of Compton camera hit energy
// ****************************************
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// yeah, yeah, I know it's ugly: hardcoded numbers (in a general library even!!!), but anyway...
namespace VIPSemiConductorConstants
{
	const double 		k_fromSigma2FWHM = 2.355;

	// constants equal for CdTe and Si
	const double 		k_Fano = 0.11;					// for CdTe and Si
	const double 		k_fit_pNoise = 3.2e-3; 			// MeV 
	const double 		k_fit_sigmaPNoise = k_fit_pNoise/k_fromSigma2FWHM;	
	
	// constants for Si
	const double 		k_epsilon_Si = 3.62e-6; 			// MeV (pair creation energy)
	
	// constants for CdTe 
	const double 		k_epsilon_CdTe = 4.43e-6; 			// MeV (pair creation energy)
	const double 		k_fit_pTrapping_CdTe = 0.014; 		// delE/E ~ 1.4% 
	const double 		k_fit_sigmaPTrapping_CdTe = k_fit_pTrapping_CdTe/k_fromSigma2FWHM;
}

// =======================================

class VIPSemiConductor
{
public: 
						VIPSemiConductor();
	virtual 			~VIPSemiConductor() {};

	virtual double 		GetSmearedEnergy( const double& in_E ) const = 0;
	virtual double 		StatisticalEnergySmearingSigma( const double& in_E ) const = 0;

protected: 

private:
	static bool			sVIPSemiConductor_RandomSeed_Done;
};

// =======================================

class VIPSi : public VIPSemiConductor
{
public:
						VIPSi() {};
	virtual 			~VIPSi() {};

	virtual double 		GetSmearedEnergy( const double& in_E ) const;
	virtual double 		StatisticalEnergySmearingSigma( const double& in_E ) const;

private:
};

// =======================================

class VIPCdTe : public VIPSemiConductor
{
public: 
						VIPCdTe() {};
	virtual 			~VIPCdTe() {};

	virtual double 		GetSmearedEnergy( const double& in_E ) const;
	virtual double 		StatisticalEnergySmearingSigma( const double& in_E ) const;

private:
};



#endif


