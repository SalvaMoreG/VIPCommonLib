
#include "CVIPPhysicsUtils.h"
#include "VIPconstants.h"

#include <cassert>
#include <cstdlib>

using namespace std;

CPhysicsUtils::CPhysicsUtils()
	: m_maxP6(0.0)
	, m_maxP3(0.0)
{
	LookUpTable( true );
	LookUpTable( false );
}

CPhysicsUtils::~CPhysicsUtils()
{
}


// *******************************************************************************************************

double
CPhysicsUtils::CalculateKleinNishinaProb( const double& in_Etot, const double& in_costh )
{
	double term = 1.0 + (in_Etot/mass_electron_keV) * (1.0 - in_costh);
	if (term == 0)
		return 0.0;

	double pratio = 1./term;

	const double k = 1.0;
	term = (pratio + 1./pratio - 1 + in_costh*in_costh)/2.0;
	double result = k * pratio * pratio * term;

	return result;
}

// *******************************************************************************************************

// *******************************************************************************************************

double
CPhysicsUtils::GetP6AbsorberPhotoElectricProbability( const double& in_E, const double& in_pathlength ) const
{
	int index; 
	double prevE, nextE;
	InterpolateEnergies( true, in_E, index, prevE, nextE );
	assert (index > 0);
	assert (int (1000 * prevE) != int (1000 * nextE) );

	// logaritmic extrapolation 
	double att1 = m_CdTePhotoElAtt[index-1];
	double att2 = m_CdTePhotoElAtt[index];
	double logcoeff = (std::log10(att2) - std::log10(att1))/(std::log10(nextE) - std::log10(prevE));
	double logpeAtt = std::log10(att1) + (std::log10(in_E) - std::log10(prevE)) * logcoeff;
	double peAtt = exp(logpeAtt * log(10));

	att1 = m_CdTeTotWiAtt[index-1];
	att2 = m_CdTeTotWiAtt[index];
	logcoeff = (std::log10(att2) - std::log10(att1))/(std::log10(nextE) - std::log10(prevE));
	double logtotAtt = std::log10(att1) + (std::log10(in_E) - std::log10(prevE)) * logcoeff;
	double totAtt = exp(logtotAtt * log(10));

	// exponential term
	// double expterm = 1.0 - exp( -1.0 * totAtt * in_pathlength );
	// double P6 = (peAtt / totAtt) * expterm;

	double expterm = exp( -1.0 * totAtt * in_pathlength );
	double P6 = peAtt * expterm;
	P6 = P6 / m_maxP6;

// assert (P6 < 1.0);
	if (P6 > 1.0) P6 = 1;

	return P6;
}

// *******************************************************************************************************

double
CPhysicsUtils::GetP3ScattererComptonProbability( const double& in_E, const double& in_pathlength ) const
{
	int index; 
	double prevE, nextE;
	InterpolateEnergies( false, in_E, index, prevE, nextE );
	assert (index > 0);
	assert (int (1000 * prevE) != int (1000 * nextE) );

	// log extrapolation 
	double att1 = m_SiComptonAtt[index-1];
	double att2 = m_SiComptonAtt[index];
	double logcoeff = (std::log10(att2) - std::log10(att1))/(std::log10(nextE) - std::log10(prevE));
	double logcompAtt = std::log10(att1) + (std::log10(in_E) - std::log10(prevE)) * logcoeff;
	double compAtt = exp(logcompAtt * log(10));

	att1 = m_SiTotWiAtt[index-1];
	att2 = m_SiTotWiAtt[index];
	logcoeff = (std::log10(att2) - std::log10(att1))/(std::log10(nextE) - std::log10(prevE));
	double logtotAtt = std::log10(att1) + (std::log10(in_E) - std::log10(prevE)) * logcoeff;
	double totAtt = exp(logtotAtt * log(10));

	// exponential term
//	double expterm = 1.0 - exp( -1.0 * totAtt * in_pathlength );
//	double P3 = (compAtt / totAtt) * expterm;

	double expterm = exp( -1.0 * totAtt * in_pathlength );
	double P3 = compAtt * expterm;

if (P3 > m_maxP3)
{
	cout << "P3 oops: " << in_E << " P3: " << P3 << " P3max: " << m_maxP3 << endl;
	cout << "prevE: " << prevE << " nextE: " << nextE << endl;
	cout << "att1: " << att1 << " att2: " << att2 << endl;
}

	P3 = P3 / m_maxP3;

assert (P3 < 1.0);

	return P3;
}

// *******************************************************************************************************


double
CPhysicsUtils::GetP5ScattererEscapeProbability( const double& in_E, const double& in_escapeLength ) const
{
	double totAtt = GetTotalAttenuation( in_E, false );

	double expterm = exp( -1.0 * totAtt * in_escapeLength );
	double P5 = expterm;
assert (P5 < 1.0);
	return P5;
}

// *******************************************************************************************************

double
CPhysicsUtils::GetTotalAttenuation( const double& in_energy, bool isCdTe ) const
{
	int index; 
	double prevE, nextE;
	InterpolateEnergies( false, in_energy, index, prevE, nextE );
	assert (index > 0);
	assert (int (1000 * prevE) != int (1000 * nextE) );

	// log extrapolation 
	double att1 = isCdTe ? m_CdTeTotWiAtt[index-1] : m_SiTotWiAtt[index-1];
	double att2 = isCdTe ? m_CdTeTotWiAtt[index] : m_SiTotWiAtt[index];

	double logcoeff = (std::log10(att2) - std::log10(att1))/(std::log10(nextE) - std::log10(prevE));
	double logtotAtt = std::log10(att1) + (std::log10(in_energy) - std::log10(prevE)) * logcoeff;
	double totAtt = std::exp(logtotAtt * log(10));

	return totAtt;
}


// *******************************************************************************************************

void
CPhysicsUtils::InterpolateEnergies( bool isCdTe, const double& in_E
								  , int& io_index, double& io_prevE, double& io_nextE ) const
{
	const std::vector<double>* energies (NULL);
	if (!isCdTe)
	{
		energies = &(m_SiEnergies);
	}
	else
	{
		energies = &(m_CdTeEnergies);
	}
	assert (energies != NULL);

    assert (energies->size() > 1);

	io_prevE = (*energies)[0];
	io_nextE = io_prevE;
	io_index = -1;
	for (int i = 1; i < energies->size() && io_index == -1 ; i++)
	{
		io_nextE = (*energies)[i];
		if (io_prevE <= in_E && in_E <= io_nextE)
			io_index = i;
		else
			io_prevE = io_nextE;
	}
	// if not found...
	if (io_index < 0)
	{
		std::string dump = (!isCdTe) ? "silicon" : " CdTe"; 
		cout << "material: " << dump << endl;
		cout << "energies: " << endl;
		cout << endl;
		cout << "prve and next: " << io_prevE << " " << io_nextE << endl;
		cout << "in_E: " << in_E << endl;
	}
	assert (io_index > 0);
	assert (int (1000 * io_prevE) != int (1000 * io_nextE) );
}

// *******************************************************************************************************

void
CPhysicsUtils::LookUpTable( bool isCdTe)
{
	std::ifstream tabfile;
	if (isCdTe)
		tabfile.open( "PHOTON_XSECTION_CdTe.dat" );
	else
		tabfile.open( "PHOTON_XSECTION_Si.dat" );

	if ( !tabfile.is_open())
	{
		cout << "ERROR. File not open: 'PHOTON_XSECTION_xxxx.dat' " << endl;
		exit(1);
	}
	// two lines of text...
	std::string dummy;
	for (int j = 0; j < 2; j++)
	{
		for (int i = 0; i < 6; i++)
		{
			tabfile >> dummy ;
		}
	}

	// here we go
	double enerMev, att_coh, att_inc, att_phot, att_totwi, att_totwo;
	double ener;

	if (isCdTe)
	{
		m_CdTeEnergies.clear();
		m_CdTePhotoElAtt.clear();
		m_CdTeTotWiAtt.clear();
	}
	else
	{
		m_SiEnergies.clear();
		m_SiComptonAtt.clear();
		m_SiTotWiAtt.clear();
	}

	// att is in 1/cm, our space units are in mm. 
	// converge att to 1/mm by multiplying with 0.1 
	double overcmtoovermm = 0.1;
	// from mass-attenuation to linear attenuation: mu_l = density * mu_mass
	double density = isCdTe ? 5.85 : 2.33;
	double factor = density*overcmtoovermm;		
	while ( !tabfile.eof() )
	{
		tabfile >> enerMev >> att_coh >> att_inc >> att_phot >> att_totwi >> att_totwo;
		if ( !tabfile.eof() )
		{
			ener = enerMev*1E3;
			att_coh = factor * att_coh;
			att_inc = factor * att_inc;
			att_phot = factor * att_phot;
			att_totwi = factor * att_totwi;
			att_totwo = factor * att_totwo;

			// all attenuations in the table are given in [cm-1]. However, we work with units in mm.
			if (isCdTe)
			{
				m_CdTeEnergies.push_back( ener );
				m_CdTePhotoElAtt.push_back( att_phot );
				m_CdTeTotWiAtt.push_back( att_totwi );

				double pathlength = 2.0;
				// double expterm = 1.0 - exp( -1.0 * att_totwi * pathlength );
				double expterm = exp( -1.0 * att_totwi * pathlength );

				// double P6 = (att_phot);		// max is als expterm == 1
				double P6 = att_phot * expterm;

				// TODO: do energy cut by interpolation....
				if ( P6 > m_maxP6 && ener > 149.0 )  
				{
					// ugly way to make sure max is always max
					m_maxP6 = P6; 
// cout << "max, ener: " << ener << " maxP6: " << m_maxP6 << endl;
				}
			}
			else
			{
				m_SiEnergies.push_back( ener );
				m_SiComptonAtt.push_back( att_inc );
				m_SiTotWiAtt.push_back( att_totwi );

				// double P3 = (att_inc/att_totwi) * expterm;
				double P3 = att_inc; // * expterm;
				if ( P3 > m_maxP3 && ener >= 499 )   
				{
					m_maxP3 = P3; 
				}
			}
		}
	}
}

// ********************************************************************************************






