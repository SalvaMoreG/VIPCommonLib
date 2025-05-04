
// IEVENT: 0
// fovidx: 12880: [0,80,0], with planeAxis: 1
// Additional: 12236: [0,76,0]  12397: [0,77,0]  12558: [0,78,0]  12719: [0,79,0]  13041: [0,81,0]  13202: [0,82,0]  13363: [0,83,0]  13524: [0,84,0]  ENDadditional

void
ParseCoordStr( const std::string& in_string, int& out_ix, int& out_iy, int& out_iz );

int 
Main()
{
	string filename;
	cout << "filename" << endl;
	cin >> filename;
	std::ifstream ffile( filename.c_str() );
	if ( !ffile.is_open() )
	{
		cout << "FILE NOT FOUND: " << filename << endl;
		return 1;
	}

	// TH3D* htest = new TH3D("htest", "TEST", 161, -20.125, 20.125, 161, -20.125, 20.125, 1, -12.5, 12.5 );
	// TH3D* htest = new TH3D("htest", "TEST", 161, -0.5, 160.5, 161, -0.5, 160.5, 1, -0.5, 0.5 );

	TH2D* hiev0 = new TH2D("hiev0", "TEST", 161, -0.5, 160.5, 161, -0.5, 160.5 );
	TH2D* hiev1 = new TH2D("hiev1", "TEST", 161, -0.5, 160.5, 161, -0.5, 160.5 );
	TH2D* hiev2 = new TH2D("hiev2", "TEST", 161, -0.5, 160.5, 161, -0.5, 160.5 );
	
	string dumStr, coordStr;
	int dumInt;
	int ix, iy, iz;
	int ievent = 0;
	while ( !ffile.eof() )
	{
		ffile >> dumStr;
		if ( !ffile.eof() )
		{
			if (dumStr == "IEVENT:")
				ffile >> ievent;
			else if (dumStr == "fovidx:")
			{
				ffile >> dumStr >> coordStr >> dumStr >> dumStr >> dumInt;
				ParseCoordStr( coordStr, ix, iy, iz );
				// htest->Fill(ix, iy, iz);
				if (ievent == 0)
					hiev0->Fill(ix, iy);
				else if (ievent == 1)
					hiev1->Fill(ix, iy);
				else if (ievent == 2)
					hiev2->Fill(ix, iy);
			}
			else if (dumStr == "Additional:")
			{
				bool domore = true;
				while (domore)
				{ 
					ffile >> dumStr;
					if (dumStr == "ENDadditional")
						domore = false;
					else
					{
						ffile >> coordStr;
						ParseCoordStr( coordStr, ix, iy, iz );
						// htest->Fill(ix, iy, iz);
						if (ievent == 0)
							hiev0->Fill(ix, iy);
						else if (ievent == 1)
							hiev1->Fill(ix, iy);
						else if (ievent == 2)
							hiev2->Fill(ix, iy);
					}
				}
			}
		}
	}

	TCanvas* m_c = new TCanvas("m_c", "TEST", 700, 700);
	m_c->Divide(2, 2);

	/*
	htest->SetMarkerColor( kRed );
	htest->SetMarkerStyle(20);
	htest->Draw("P");
	*/
	m_c->cd(1);
	hiev0->SetMarkerColor( kRed );
	hiev0->SetMarkerSize( 0.5 );
	hiev0->SetMarkerStyle(20);
	hiev0->Draw("P");
	m_c->cd(2);
	hiev1->SetMarkerColor( kRed );
	hiev1->SetMarkerSize( 0.5 );
	hiev1->SetMarkerStyle(20);
	hiev1->Draw("P");
	m_c->cd(3);
	hiev2->SetMarkerColor( kRed );
	hiev2->SetMarkerSize( 0.5 );
	hiev2->SetMarkerStyle(20);
	hiev2->Draw("P");

}

// [0,76,0]  
void
ParseCoordStr( const std::string& in_string, int& out_ix, int& out_iy, int& out_iz )
{
	int posx1 = in_string.find(string("["));
	int posx2 = in_string.find(string(","));
	int len = posx2 - posx1 - 1;
	string xstr = in_string.substr ( posx1+1, len);
	out_ix = atoi(xstr.c_str());

	len = in_string.size();
	string remains = in_string.substr( posx2+1, len);
	int posy1 = 0;
	int posy2 = remains.find(string(","));
	len = posy2 - posy1;
	string ystr = remains.substr ( posy1, len);
	out_iy = atoi(ystr.c_str());

	len = remains.size();
	remains = remains.substr( posy2+1, len);
	int posz1 = 0;
	int posz2 = remains.find(string("]"));
	len = posz2 - posz1;
	string zstr = remains.substr ( posz1, len);
	out_iz = atoi(zstr.c_str());

	// cout << "in_string: " << in_string << " x: " << xstr << " y: " << ystr << " z: " << zstr << endl;
}


