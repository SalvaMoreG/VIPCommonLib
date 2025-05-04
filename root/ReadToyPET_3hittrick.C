
int Main()
{
	string fname;
	cout << "give filename" << endl;
	cin >> fname;
	ifstream ffile( fname.c_str() );
	if ( !ffile.is_open() )
	{
		cout << "file does not exist: " << fname << endl;
		return 1;
	}

    double hit1_x, hit1_y, hit1_z, time1 ;
    double hit2_x, hit2_y, hit2_z, time2 ;
    double hit3_x, hit3_y, hit3_z, time3  ;
	double sourcePos_x, sourcePos_y, sourcePos_z; 
	double sourceGuess_x, sourceGuess_y, sourceGuess_z; 
    double srcError;

    TH1D* h1 = new TH1D("h1", "Distance S - S' (mm)", 100, -50., 50.);
	
    int alternate(-1);
	while ( !ffile.eof() )
	{
        ffile >> hit1_x >> hit1_y >> hit1_z >> time1
              >> hit2_x >> hit2_y >> hit2_z >> time2
              >> hit3_x >> hit3_y >> hit3_z >> time3
              >> sourcePos_x >> sourcePos_y >> sourcePos_z
              >> sourceGuess_x >> sourceGuess_y >> sourceGuess_z
			  >> srcError;

		if ( !ffile.eof() )
		{
            h1->Fill(srcError*alternate);   // Get a Gaussian shape
            alternate *= -1;
		}
	}

	h1->Draw();

	return 0;
}

//
/*

- Read file
- Calculate S'
    *) center = (hit_A + hit_B)/2
    *) delta_t_AB = (t_A - t_B)
    *) Delta = c * delta_t_AB/2
    *) direction = (hit_A - hit_B)/|(hit_A - hit_B)|
    *) SO: guessedSrc = center - direction * Delta
- Other stuff:
    *) c = A_S/dt_A = F_S/dt_F = B_S/dt_B
       So:
          dt_A = A_S/c
          dt_B = B_S/c
          dt_F = F_S/c
    *) delta_t_AF =  (t_A - t_F)
    *) delta_t_BF =  (t_B - t_F)
- What we know (from measurement):
	- A_B
	- A_S', B_S' and F_S'
	- (t_A - t_B)'
	- (t_A - t_F)'
	- (t_B - t_F)'
- What we can calculate (c = speed of light):
	- dt_F' = (F_S')/c
	- 
*/
