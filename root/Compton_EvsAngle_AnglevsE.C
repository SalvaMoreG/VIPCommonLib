bool
GetComptonAngleDegrees(const double& E1_keV, const double& Etot_keV, double& out_angleComptonDegrees );

bool
GetComptonEnergyFromAngleRad(const double& in_angleComptonRadians, const double& Etot_keV, double& out_E1_keV );

const double kPI = 3.14159265359;
const double mass_electron_keV = 511.0;

int Main()
{
	cout << "Give total gamma energy (keV) (511, 1274, ...)" << endl;
	double Etot_keV; cin >> Etot_keV;

	TH1D* h1 = new TH1D("h1", "Energy vs Angle", 730, 0., 185.0 );
	TH1D* h2 = new TH1D("h2", "Angle vs Energy", 2*Etot_keV, 0., Etot_keV );

	for (int ia = 0; ia < 1800; ia++)
	{
		double angleDeg = 0.1 * (double) (ia);
		double angleRad = kPI * angleDeg/180.0;
		double E1_keV;
		bool ok = GetComptonEnergyFromAngleRad( angleRad, Etot_keV, E1_keV );
		if (ok)
		{
			h1->Fill( angleDeg, E1_keV);

			if ( fabs(angleDeg - 180.0) < 0.15 )
			{
				cout << "Compton edge, angle: " << angleDeg << " degrees, E1: " << E1_keV << endl;
			}
		}
	}

	for (int iE = 0; iE < 2*Etot_keV; iE++)
	{
		double E1_keV = 0.5 * (double) (iE);
		double angleDeg;
		bool ok = GetComptonAngleDegrees( E1_keV, Etot_keV, angleDeg );
		if (ok)
		{
			h2->Fill( E1_keV, angleDeg);
		}
	}

	TCanvas* m_c1 = new TCanvas("m_c1", " ", 1200, 600);
	m_c1->Divide(2, 1);
	m_c1->cd(1);
	h2->Draw();
	m_c1->cd(2);
	h1->Draw();
}

//	=========================================================================================

bool
GetComptonAngleDegrees(const double& E1_keV, const double& Etot_keV, double& out_angleComptonDegrees )
{
   out_angleComptonDegrees = -1000;

   if (Etot_keV > 0 && (Etot_keV - E1_keV) > 0)
   {
       double coscompton = 1.0 - mass_electron_keV*(1.0/(Etot_keV - E1_keV) - 1.0/Etot_keV);
       if (fabs(coscompton) <= 1)
       {
          out_angleComptonDegrees = acos(coscompton) * 180.0/kPI;
          return true;
       }
   }
   return false;
}

bool
GetComptonEnergyFromAngleRad(const double& in_angleComptonRadians, const double& Etot_keV, double& out_E1_keV )
{
	double cosTh = cos(in_angleComptonRadians);
	double tmp = 1.0 + mass_electron_keV / (Etot_keV * (1.0 - cosTh)) ;

	out_E1_keV = Etot_keV / tmp;
	
	// Since all angles are permitted (0 - 180 degrees; 
	//	higher than 180 degrees gives same cosine as smaller than 180 degrees)
	// this function always returns "true"
	return true;
}

// ======================================================================

