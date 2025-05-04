#include <iostream>
#include <fstream>

#include <string>
#include <cmath>

#include "CVIPRandom.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

using namespace std;

//  int algorithm_poisson_random_number(const double& in_lambda);

//
// https://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
//
double GetNextTime(const double& in_lambda, double& io_random);
double GetRandom();

int main()
{
	const int ntrials(100000);

	// double lambda(1.0); 	// 1 Hz

	double lambda(2341.0);	 // 200 Hz

	// double lambda(200.0e-6); // 200*10^-6 events per 1 us

	// double lambda(200.0e-9); // 200*10^-9 events per 1 ns;
    
    cout << "give lambda (if 0, use old value): " << lambda << endl;
    double wtv(0.0) ; cin >> wtv;
    if (wtv > 0.0)
        lambda = wtv;
    
    cout << "FINAL lamba: " << lambda << endl;

	TTree* timeTree = new TTree("timedifferences", "TIMEDIFFERENCES");
	double m_deltat; 
	timeTree->Branch("deltat", &m_deltat);

	TH1D* hr = new TH1D("hr", "RANDOMS", 1001, 0., 1.01);

	string dummy17("EBBD_42_34000");
	string dummy8("AA67_42_34000");
	double energy_MeV(0.511), x(60.0), y(60.0), z(0.0);

	ofstream outfile("LM_TIME_NOT_ORDERED.dat_NEW");

	unsigned long long timestamp1(0);
	unsigned long long timestamp2(0);

	int icount = 0;
	for (int itrial = 0; itrial < ntrials; itrial++)
	{
        // One chip
		for (int iloop = 0; iloop < ntrials/100; iloop++)
		{
			double random;
			m_deltat = GetNextTime(lambda, random);
			timeTree->Fill();

			timestamp1 += m_deltat;
			hr->Fill(random);
			outfile << itrial << " " << dummy17 << " " << timestamp1*1000 << " "
					<< energy_MeV << " " << x << " " << y << " " << z << endl;
		}

		// Another chip
		for (int iloop = 0; iloop < ntrials/100; iloop++)
		{
			double random;
			m_deltat = GetNextTime(lambda, random);
			timeTree->Fill();

			timestamp2 += m_deltat;
			hr->Fill(random);
			outfile << itrial << " " << dummy8 << " " << timestamp2*1000 << " "
					<< energy_MeV << " " << x << " " << y << " " << z << endl;

			itrial++;
		}

		timestamp1 += 10e9;     // skip 10e9 ns ( = 1 second)
		timestamp2 += 10e9;
	}

	TFile aFile("poisson.root", "RECREATE");

	timeTree->Write();
	hr->Write();

	return 0;
}

// ===============================================================================
/*
int algorithm_poisson_random_number(const double& in_lambda)
{
    // init: //
	double L = std::exp(-1. * in_lambda);
	int k = 0;
	double p = 1.0;
	while (p > L)
	{
    	k = k + 1;
        // Generate uniform random number u in [0,1] and let p ← p × u.
		double u = GetRandom(); // rand() / (RAND_MAX + 1.);
		p = p * u;
	}
    // while p > L.
	return k - 1;
}
*/
// ===============================================================================

// See: https://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
//
// 	def nextTime(rateParameter):
// 	    return -math.log(1.0 - random.random()) / rateParameter

double GetNextTime(const double& in_lambda, double& io_random)
{
	io_random = GetRandom();    // uniform value between 0 and 1
	double nexttime = -1.0 * std::log(1.0 - io_random) / in_lambda;
    // cout << io_random << " " << nexttime << endl;
	return nexttime;
}


// ===============================================================================

double GetRandom()
{
	/*
	bool isok(false);
	while (!isok)
	{
		io_random = rand() / (RAND_MAX + 1.0);
		isok = (io_random >= 0 && io_random < 1.0);
	}
	*/

	CVIPRandomUniform uniformDist(0.0, 1.0);
	double X = uniformDist.GetNewValue();

	return X;
}




