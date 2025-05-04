
#include "VIPconstants.h"
// 
#include <random>
#include <iostream>
#include <fstream>

#include "TH1D.h"
#include "TCanvas.h"

double getangle(const double& in_x, const double& in_y);

void FillCylinderPETXY( const double& in_delR, const double& in_delT, 
    TH1D* in_hR, TH1D* in_hT, TH1D* in_hX, TH1D* in_hY,
    double& io_x, double& io_y);

using namespace std;

// ================================================================

int main(int argc, char* argv[])
{
	string arg;
	bool onlyZ(false);
	bool onlyXY(false);
	bool doGauss(false);
	int argidx(1);
    double deltaZ(1.0);
    double deltaXY(0.5);
    bool planarXY(false);
	if (argc > 1)
	{
		for (int argidx = 1; argidx < argc; argidx++)
		{
			arg = argv[argidx];
			if (arg == "-onlyZ")
			{
				onlyZ = true;
				cout << "Only randomizing in Z direction" << endl;
			}
			else if (arg == "-deltaZ")
			{
                if ( (argidx+1) >= argc)
                {
                    cout << "Wrong number of arguments. Usage: " << endl;
                    cout << argv[0] << " -deltaZ <deltaZ> " << endl;
                    exit(1);
                }
                argidx++;
                deltaZ = atof( argv[argidx] );
			}
			else if (arg == "-deltaXY")
			{
                if ( (argidx+1) >= argc)
                {
                    cout << "Wrong number of arguments. Usage: " << endl;
                    cout << argv[0] << " -deltaXY <deltaXY> " << endl;
                    exit(1);
                }
                argidx++;
                deltaXY = atof( argv[argidx] );
			}
			else if (arg == "-planarXY")
			{
				planarXY = true;
				cout << "PET XY planes (no ring)" << endl;
			}
			else if (arg == "-onlyXY")
			{
				onlyXY = true;
				cout << "Only randomizing in XY direction" << endl;
			}
			else if (arg == "-gauss")
			{
				doGauss = true;
				cout << "using normal distribution (Gauss)" << endl;
			}
			else 
			{
				cout << argv[0] << ": randomize in X, Y[+-0.5mm] and Z[+-1mm]" << endl;
                cout << argv[0] << " -deltaZ <deltaZ> " << endl;
                cout << argv[0] << " -deltaXY <deltaXY> " << endl;
				cout << argv[0] << " -onlyZ: randomize only in Z[+-1mm]" << endl;
                cout << argv[0] << " -planarXY: XY are not along a ring (PET = two planes)" << endl;
				cout << argv[0] << " -onlyXY: randomize only in X and Y[+-0.5mm]" << endl;
				cout << argv[0] << " -gauss: using normal distribution (Gauss)" << endl;
				return 0;
			}
		}
	}

    if ( !planarXY )
	{
		cout << "WARNING!!! It's better to run with flag -planarXY!!!" << endl;
	}

	std::string fname;
	cout << "filename?" << endl;
	cin >> fname;

	std::ifstream ffile( fname.c_str() );
	if (!ffile.is_open())
	{
		cout << "File not open: " << fname << endl;
		return 1;
	}

	// UNIFORM distribution
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> distXY(-1.0*deltaXY, deltaXY);
    std::uniform_real_distribution<double> distZ(-1.0*deltaZ, deltaZ);
    // X = dist(mt);

	// NORMAL distribution
    std::random_device rd2;
    std::mt19937 mt2(rd2());
	double sigmaXY = deltaXY/3;	// 3 sigma = 99%
    std::normal_distribution<double> distGaussXY(0.0, sigmaXY);
	double sigmaZ = deltaZ/3;	// 3 sigma = 99%
    std::normal_distribution<double> distGaussZ(0.0, sigmaZ);
    // X = distribution( mt );

	ofstream outfile("out_COINC_PETLORS_RANDOMIZED.dat_NEW");

    TH1D* hR = new TH1D("hR", "delta R distribution check", 1000, -2.0*deltaXY, 2.0*deltaXY);
    TH1D* hT = new TH1D("hT", "delta T distribution check", 1000, -2.0*deltaXY, 2.0*deltaXY);    
    TH1D* hX = new TH1D("hX", "delta X distribution check", 1000, -2.0*deltaXY, 2.0*deltaXY);
    TH1D* hY = new TH1D("hY", "delta Y distribution check", 1000, -2.0*deltaXY, 2.0*deltaXY);
	TH1D* hZ = new TH1D("hZ", "delta Z distribution check", 1000, -1.5*deltaZ, 1.5*deltaZ);

	double z1, y1, x1, e1, z2, y2, x2, e2;
    double delR, delT, delZ;
	while ( !ffile.eof() )
	{
		ffile >> z1 >> y1 >> x1 >> e1 >> z2 >> y2 >> x2 >> e2;
		ffile.ignore(1024, '\n');

		if ( !ffile.eof() )
		{
			// GAMMA 1 --------------------------------------------
			if (!onlyZ)  // XY 1
			{
				delR = (doGauss) ? distGaussXY(mt2) : distXY(mt);
				delT = (doGauss) ? distGaussXY(mt2) : distXY(mt);
                
                if (planarXY)		// simple X and Y variation
                {
                    x1 += delR;
                    y1 += delT;
                    hX->Fill(delR);
                    hY->Fill(delT);
                }
                else				// vary X and Y with radial and tangential components
                    FillCylinderPETXY( delR, delT, hR, hT, hX, hY, x1, y1 );
                
			}
			if (!onlyXY)  // Z 1
			{
				delZ = (doGauss) ? distGaussZ(mt2) : distZ(mt);
				z1 += delZ;
				hZ->Fill(delZ);
			}

			// GAMMA 2 --------------------------------------------
			if (!onlyZ)  // XY 2
			{
				delR = (doGauss) ? distGaussXY(mt2) : distXY(mt);
				delT = (doGauss) ? distGaussXY(mt2) : distXY(mt);
                
                if (planarXY)		// Simple X and Y variation
                {
                    x2 += delR;
                    y2 += delT;
                    hX->Fill(delR);
                    hY->Fill(delT);
                }
                else				// vary X and Y with radial and tangential components
                    FillCylinderPETXY( delR, delT, hR, hT, hX, hY, x2, y2 );

			}
			if (!onlyXY)  // Z2
			{
            	delZ = (doGauss) ? distGaussZ(mt2) : distZ(mt);
				z2 += delZ;
				hZ->Fill(delZ);
			}

			outfile << z1 << " " << y1 << " " << x1 << " " << e1 << " " 
                    << z2 << " " << y2 << " " << x2 << " " << e2 << endl;
		}
	}

	TCanvas* m_1 = new TCanvas("m_1", " ", 2400, 1600);
    m_1->Divide(3, 2);
    m_1->cd(1);
    hR->Draw("HIST");
    m_1->cd(2);
    hT->Draw("HIST");    
    m_1->cd(3);
	hZ->Draw("HIST");    
    m_1->cd(4);
    hX->Draw("HIST");
    m_1->cd(5);
    hY->Draw("HIST");    
	m_1->Print("checkRandom.png", "png");
}

// ========================================================

double getangle(const double& in_x, const double& in_y)
{
	double tanalpha, alpharad;
	int quadrant(0);
	if (in_x >= 0 && in_y > 0) quadrant = 1;
	else if (in_x < 0 && in_y > 0) quadrant = 2;
	else if (in_x <= 0 && in_y < 0) quadrant = 3;
	else if (in_x > 0 && in_y < 0) quadrant = 4;
	if (in_x != 0)
	{
		tanalpha = in_y/in_x;
		alpharad = atan(tanalpha);
		if (quadrant == 2 || quadrant == 3) alpharad += kPI;
	}
	else
	{
		alpharad = (quadrant == 1) ? 0.5 * kPI : -0.5 * kPI;
	}
	return alpharad;
}

// ========================================================

void FillCylinderPETXY( const double& in_delR, const double& in_delT, 
    TH1D* in_hR, TH1D* in_hT, TH1D* in_hX, TH1D* in_hY,
    double& io_x, double& io_y)
{
    double delX_r, delY_r, delX_t, delY_t, anglerad;
    
    anglerad = getangle(io_x, io_y); 
    delX_r = in_delR * cos(anglerad);
    delY_r = in_delR * sin(anglerad);
    delY_t = in_delT * cos(anglerad);
    delX_t = in_delT * sin(anglerad);
    io_x += delX_r + delX_t;
    io_y += delY_r + delY_t;
    in_hR->Fill(in_delR);
    in_hT->Fill(in_delT);
    in_hX->Fill(delX_r + delX_t);
    in_hY->Fill(delY_r + delY_t);
}


