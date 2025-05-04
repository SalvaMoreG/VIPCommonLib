
double CalculateDistance( const double& in_x1, const double& in_y1
						, const double& in_x2, const double& in_y2 );

const double k_cspeed(3E-1); // mm per picoseconds

int Main()
{
	double A_x = -20.0;
	double A_y =   1.0;
	double B_x =  18.5;
	double B_y =   6.3;
	double F_x =  19.0;
	double F_y =  -4.0;

	double delta_T_FA_clean = -56.32;	
	double delta_X_FA_clean = delta_T_FA_clean * k_cspeed;
	double delta_T_FB_clean =  14.45;
	double delta_X_FB_clean = delta_T_FB_clean * k_cspeed;

	double len_AB_lor = CalculateDistance( A_x, A_y, B_x, B_y );
	double len_AF = CalculateDistance( A_x, A_y, F_x, F_y );
	double len_BF = CalculateDistance( B_x, B_y, F_x, F_y );
	cout << "LOR: " << len_AB_lor << endl;
	
	double lor_from_B_to_A_x = A_x - B_x;
	double lor_from_B_to_A_y = A_y - B_y;
	double lor_B2A_unit_x = lor_from_B_to_A_x/len_AB_lor;
	double lor_B2A_unit_y = lor_from_B_to_A_y/len_AB_lor;
	cout << "LOR unit: " << lor_B2A_unit_x << " , " << lor_B2A_unit_y << endl;

	// Going from A to B
	double AS = 1;
	double error(0);
	double preverror = 10000;
	double prev_S_x(0), prev_S_y(0);
	bool stop(false);
	while (AS < len_AB_lor && !stop)
	{
		double S_x = A_x + (-1. * AS * lor_B2A_unit_x);
		double S_y = A_y + (-1. * AS * lor_B2A_unit_y);

		double len_FS = CalculateDistance( S_x, S_y, F_x, F_y );
		double delta_T_AS = AS / k_cspeed;
		double delta_T_FS = len_FS / k_cspeed;
		double guess_delta_T_FA = delta_T_FS - delta_T_AS; 
		
		double error = fabs(guess_delta_T_FA - delta_T_FA_clean);
		if (error < 2)
		{
			cout << "AS: " << AS << " FS: " << len_FS << " guess T_FA: " << guess_delta_T_FA 
			     << " S: " << S_x << " " << S_y << " error: " << error << endl;
		}
		if (error < preverror)
		{
			prev_S_x = S_x;
			prev_S_y = S_y;
			preverror = error;
		}
		else
			stop = true;
		AS += 0.5;
	}
	cout << "Final error: " << preverror << endl;
	cout << "final source (AF): " << prev_S_x << " , " << prev_S_y << endl;

	// Going from B to A
	double BS = 1;
	error = 0;
	preverror = 10000;
	prev_S_x = 0, prev_S_y = 0;
	stop = false;
	while (AS < len_AB_lor && !stop)
	{
		double S_x = B_x + (BS * lor_B2A_unit_x);
		double S_y = B_y + (BS * lor_B2A_unit_y);

		double len_FS = CalculateDistance( S_x, S_y, F_x, F_y );
		double delta_T_BS = BS / k_cspeed;
		double delta_T_FS = len_FS / k_cspeed;
		double guess_delta_T_FB = delta_T_FS - delta_T_BS; 
		
		double error = fabs(guess_delta_T_FB - delta_T_FB_clean);
		if (error < 5)
		{
			cout << "BS: " << BS << " FS: " << len_FS << " guess T_FB: " << guess_delta_T_FB 
			     << " S: " << S_x << " " << S_y << " error: " << error << endl;
			cout << "	" << " dT_BS: " << delta_T_BS << " dT_FS: " << delta_T_FS << endl;
		}
		if (error < preverror)
		{
			prev_S_x = S_x;
			prev_S_y = S_y;
			preverror = error;
		}
		else
			stop = true;
		BS += 0.1;
	}
	cout << "Final error: " << preverror << endl;
	cout << "final source (BF): " << prev_S_x << " , " << prev_S_y << endl;
	
}

double CalculateDistance( const double& in_x1, const double& in_y1
						, const double& in_x2, const double& in_y2 )
{
	double d_x = in_x1 - in_x2;
	double d_y = in_y1 - in_y2;
	double d = sqrt(d_x * d_x + d_y * d_y);
	return d;
}



