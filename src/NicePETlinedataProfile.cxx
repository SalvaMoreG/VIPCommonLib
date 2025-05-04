
#include <stdio.h>
#include <fstream>
#include <iostream>

#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TError.h"
#include "TLine.h"
#include "TText.h"

using namespace std;

enum ALGOTYPE
{
	ALG_NONE = 0,
	ALG_OSEM,
	ALG_FBP, 
	ALG_OE
};

ALGOTYPE m_algoType( ALG_NONE );

bool m_isFromKatya = false;

void MakeLineProfilePlot( const std::string& in_fname );
void PrintCanvas( const std::string in_fname, TCanvas* in_c );

// =====================================================

int main( int argc, char *argv[] )
{
	std::string fname;
	std::string flistname;
	cout << "give name of the file containing all tsv files" << endl;
	cin >> flistname;

	std::ifstream flist( flistname.c_str() );
	while ( !flist.eof() )
	{
		flist >> fname;
		if ( !flist.eof() )	
		{
			cout << "dealing with file: " << fname << endl;
			MakeLineProfilePlot( fname );
		}
	}

	return 0;
}

// =====================================================

void
MakeLineProfilePlot( const std::string& in_fname )
{
	double k,l;
	ifstream fdata;

	fdata.open( in_fname.c_str(), ios::in);
	if (!fdata.is_open()) 
	{
		cout << "ERROR: File not open: " << in_fname << endl;
	}
	cout << "Plotting file: " << in_fname << endl;

	// oh jeez....
	std::string aLine;

	bool isHeader(true);
	const string searchStr("#");
	streampos oldpos( 0 );
	int numHlines( 0 );
	while (isHeader)
	{
		oldpos = fdata.tellg();  // stores the position
		std::getline ( fdata, aLine);
		size_t position = aLine.find( searchStr );
		// oldpos = fdata.tellg();  // stores the position
		if (position == string::npos)
		{
            // cout << "read line: " << aLine << endl;
			isHeader = false;
			fdata.seekg (oldpos);
		}
		else
			numHlines++;
	}
	cout << "Read " << numHlines << " comment lines." << endl;

	int numDlines = 0;
	while (!fdata.eof()) 
	{
		fdata >> k >> l;
		if (numDlines == 0) 
			cout << "first values: " << k << " " << l << endl;
		if (fdata.eof()) 
			cout << "last values: " << k << " " << l << endl;

		if (!fdata.eof())
			numDlines++;
 	}
	cout << "number of points: " << numDlines << endl;

	fdata.close();
	// opening again
	fdata.open( in_fname.c_str(),ios::in);

	// Reading header again
	for (int j = 0; j < numHlines; j++)
	{
		std::getline ( fdata, aLine);
		// cout << "line " << j << " : " << aLine << endl;
	}
	
	double* x = new double[1000];
	double* y = new double[1000];
	for (int i = 0; i < numDlines; i++)
	{
		fdata >> k >> l;
		// cout << "i: " << i << " x: " << k << " y: " << l << endl;
		x[i] = k;
		y[i] = l;
	}

	TStyle *style = new TStyle("gdlPRL","dry style for PRL");
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetFrameBorderMode(0);
	style->SetFrameFillStyle(0);
	style->SetCanvasColor(0);

	style->SetPadColor(0);
	style->SetOptStat(1);
    // style->SetOptStat(1);
	style->SetTextFont(132);

	style->SetTitleFont(132);
	style->SetTitleFontSize(0.08);
	style->SetTitleFillColor(0);
	style->SetTitleBorderSize(0);

	style->SetLabelFont(132,"X");
	style->SetLabelFont(132,"Y");

	style->SetTitleFont(132,"X");
	style->SetTitleOffset(0.7,"X");

	style->SetLabelSize(0.06,"X");
	style->SetLabelSize(0.06,"Y");

	style->SetTitleFillColor(0);
	style->SetTitleXSize(0.07);
	style->SetTitleYSize(0.07);

	gROOT->SetStyle("gdlPRL");

	TCanvas* c1 = new TCanvas("c1","c1", 600, 600);
	c1->SetGridx();
	c1->SetGridy();

	TGraph* gl4 = new TGraph(numDlines, x, y);

	if (m_algoType == ALG_NONE)
		gl4->SetTitle(" ");
	else if (m_algoType == ALG_OSEM)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. OSEM");
	else if (m_algoType == ALG_FBP)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. FBP");
	else if (m_algoType == ALG_OE)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. OE");

	gl4->SetLineColor(kBlue);
	gl4->SetLineWidth(2);
	gl4->GetXaxis()->SetTitle("Distance (mm)");

	if (m_algoType == ALG_OSEM)
		gl4->SetMaximum(0.2);	// OSEM
	else if (m_algoType == ALG_FBP)
		gl4->SetMaximum(200);	// FBP 
	else if (m_algoType == ALG_OE)
		gl4->SetMaximum(80000);	// OE 

	style->SetOptStat(1);
	gl4->DrawClone("AL");
	cout << "fit gaus? (1/0)?" << endl;
	bool dumdoit; cin >> dumdoit;
	if (dumdoit)
	{
		gl4->Draw("AL");
		gStyle->SetOptFit(1);
		gl4->Fit("gaus");
	}
	else
		style->SetOptStat(1);

	fdata.close();

    /*
	bool drawlines ( false );
	// cout << "Do you want to draw horizontal lines? (1/0)" << endl;
	// cin >> drawlines;
	if (drawlines)
	{
		TLine lineAve;
		lineAve.SetLineColor( kRed );
		lineAve.SetLineWidth( 2 );
		TText lineText;

		double value, SNR, sigh;
		int argh;

		cout << "give <peak value>" << endl;
		cin >> value;
		lineAve.SetX1(x[0]);
		lineAve.SetX2(x[numDlines-1]);
		lineAve.SetY1(value);
		lineAve.SetY2(value);
		lineAve.DrawClone();

		SNR = value;

		argh = int (10 * value + 0.5);
		sigh = argh/10.0;
		
		TString str1("<peak>: ");
		str1 += sigh;
		lineText.SetText(x[0], value, str1 );
		lineText.DrawClone();

		double SNRypos = 0.5 * value;

		cout << "give <valley>" << endl;
		cin >> value;
		lineAve.SetX1(x[0]);
		lineAve.SetX2(x[numDlines-1]);
		lineAve.SetY1(value);
		lineAve.SetY2(value);
		lineAve.DrawClone();

		SNR = SNR/value;

		argh = int (10 * value + 0.5);
		sigh = argh/10.0;
		TString str2("<valley>: ");
		str2 += sigh;
		lineText.SetText(x[0], value, str2 );
		lineText.DrawClone();

		SNR = (int (10 * SNR) + 0.5)/10.0;
		TString str3("SNR: ");
		str3 += SNR;
		lineText.SetText(x[0], SNRypos, str3);
		lineText.DrawClone();
	}
	*/

	PrintCanvas( in_fname, c1 );

	return;
}

// ==========================================================

void
PrintCanvas( const std::string in_fname, TCanvas* in_c )
{
	std::string outname;
	int pos = in_fname.find(string(".tsv"));
	int len = pos;
	if (pos != std::string::npos)
	{
		outname = in_fname.substr ( 0, len );
		outname += ".png";
	}
	else
		outname = in_fname;

	in_c->Print( outname.c_str(), "png");
}











