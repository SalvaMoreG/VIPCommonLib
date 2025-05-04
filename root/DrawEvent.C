
// e1: 512.146 x1: -22.6711 y1: -119.392 z1: -40.2805 e2: 509.564 x2: 100.255 y2: 80.4546 z2: 94.7445


void Main()
{
	gStyle->SetOptStat(0);

	double x1 = -22.6711;
	double y1 = -119.392;
	double z1 = -40.2805; 
	double x2 = 100.255; 
	double y2 = 80.4546; 
	double z2 = 94.7445;

    cout << "Give Hit 1 (x, y, z)" << endl;
    cin >> x1 >> y1 >> z1;
    cout << "Give Hit 2 (x, y, z)" << endl;
    cin >> x2 >> y2 >> z2;

	TCanvas* m1 = new TCanvas("m1","XY", 700,700);
	TH2D* h1 = new TH2D("h1", " XY frame", 250, -125, 125, 250, -125, 125); 
	h1->Draw();

    /*
	TMarker3DBox* FOVbig = new TMarker3DBox(0, 0, 0, 0.5*fovSizeX, 0.5*fovSizeY, 0.5*fovSizeZ, 0, 0 );
	FOVbig->SetLineColor(kRed);
	FOVbig->Draw("SAME");
	*/

	TLine* aLine = new TLine(x1, y1, x2, y2);
	aLine->SetLineColor(kBlue);
	aLine->SetLineWidth(3);
	aLine->Draw("SAME");
    
    TEllipse* el1 = new TEllipse(0.0,0.0, 50.0, 50.0);
    el1->SetLineColor(kRed);
    el1->SetFillStyle(0);
    el1->Draw("SAME");
    
    
    // FOV Box, change according to what you need...
    TBox* FOV = new TBox(-75, -75, 75, 75);
    FOV->SetFillStyle(0);
    FOV->SetLineColor(kGreen);
    FOV->Draw("SAME");
    


    /*
	TCanvas* m2 = new TCanvas("m2","XZ", 700,700);
	TH2D* h2 = new TH2D("h2", " XZ frame", 250, -125, 125, 250, -125, 125); 
	h2->Draw();

	TLine* aLine2 = new TLine(x1, z1, x2, z2);
	aLine2->SetLineColor(kBlue);
	aLine2->SetLineWidth(3);
	aLine2->Draw("SAME");
	*/
    

}


