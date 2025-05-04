/*
Use ConvertImgToWeighted3DTree to convert "out_binary_amide_ITER4.img_NEW" into "3D_TREE.root"
*/

int Main()
{
	cout << "Give name" << endl;
	TString fname;
	cin >> fname;

	TFile* file0 = TFile::Open(fname);
	TTree* the3DTree = (TTree*)file0->Get("the3DTree"); 

	gStyle->SetOptStat(0);

	TCanvas* m_1 = new TCanvas("m_1", " ", 2000, 1000);
	m_1->Divide(2, 1);

	m_1->cd(1);
	TH1D* h1 = new TH1D("h1", "X Caudates", 180, -90., 90.);
	the3DTree->Draw("x>>h1","WEIGHT*(z > -30 && z < -20 && y > -34 && y < -26)","HIST");
	h1->SetLineWidth(2);

	m_1->cd(2);
	TH1D* h2 = new TH1D("h2", "X Putamens", 180, -90., 90.);
	the3DTree->Draw("x>>h2","WEIGHT*(z > -30 && z < -20 && y > -4 && y < 4)","HIST");
	h2->SetLineWidth(2);

	return 0;
}

