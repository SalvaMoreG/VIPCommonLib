
#include "CVIPSemiConductorEnergy.h"

#include "CVIPRandom.h"
#include <cmath>

using namespace std;

bool VIPSemiConductor::sVIPSemiConductor_RandomSeed_Done = false;

VIPSemiConductor::VIPSemiConductor()
{
	/*
	if (!sVIPSemiConductor_RandomSeed_Done)
	{
#ifdef USING_CLHEP
		sVIPSemiConductor_RandomSeed_Done = true;
		time_t systime = time(NULL);
		long seed = (long) systime;
		CLHEP::HepRandom::setTheSeed( seed );
#endif
	}
	*/
}

// ==============================

double
VIPSi::GetSmearedEnergy( const double& in_E ) const
{
	double sigmaStat = StatisticalEnergySmearingSigma( in_E );
        // Small amount of stat. variation (# of charge carriers produced)

	CVIPRandomGauss gaussDistE(in_E, sigmaStat);
	double newE = gaussDistE.GetNewValue();
	CVIPRandomGauss gaussDistDelE(0., VIPSemiConductorConstants::k_fit_sigmaPNoise);
	double delE = gaussDistDelE.GetNewValue();

	newE = newE + delE;
       // lifetimes in Silicon are so high (compared with CdTe),
       //      that we can ignore the trapping part of the energy resolution

	return newE;
}

double
VIPSi::StatisticalEnergySmearingSigma( const double& in_E ) const
{
	double sigma = std::sqrt( in_E * VIPSemiConductorConstants::k_Fano * VIPSemiConductorConstants::k_epsilon_Si );
	return sigma;
}

// ====================================

double
VIPCdTe::GetSmearedEnergy( const double& in_E ) const
{
	// Statistical Smearing (FANO factor etc....)   --> newE
    // Electrical noise (FANO factor etc....)  -> delE
    // Trapping, tricky..... --> factor

	double sigmaStat = StatisticalEnergySmearingSigma( in_E );	
		// Small amount of stat. variation (# of charge carriers produced)

	CVIPRandomGauss gaussDistE(in_E, sigmaStat);
	double newE = gaussDistE.GetNewValue();
	CVIPRandomGauss gaussDistDelE(0., VIPSemiConductorConstants::k_fit_sigmaPNoise);
	double delE = gaussDistDelE.GetNewValue();
	CVIPRandomGauss gaussDistFactor(1., VIPSemiConductorConstants::k_fit_sigmaPTrapping_CdTe);
	double factor = gaussDistFactor.GetNewValue();

	newE = newE + delE;
    // Trapping, constant energy-factor (some part of the charge carriers are lost)

	newE = factor * newE;
       // TODO: in fact, this is not completely correct, this should be the fraction subtracted from the energy...
       // But result is more or less ok

	return newE;
}

double
VIPCdTe::StatisticalEnergySmearingSigma( const double& in_E ) const
{
	double sigma = std::sqrt( in_E * VIPSemiConductorConstants::k_Fano * VIPSemiConductorConstants::k_epsilon_CdTe );
	return sigma;
}

// =======================================


