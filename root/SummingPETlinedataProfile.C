#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

bool IsHeader(const string& in_Line);
int GetHistoParameters(const string& fname, double& xmin, double& xmax);

// ============================================

int Main()
{
	string flistname;	
	cout << "give filename of list containing all tsv files" << endl;
	cin >> flistname;
	ifstream flist( flistname.c_str(), ios::in);
	if ( !flist.is_open() )
	{
		cout << "File does not exist: " << flistname << endl;
		return 1;
	}
	string fname("");
	int nbins; 
	double xmin, xmax;
	while ( !flist.eof() )
	{
		flist >> fname;
		nbins = GetHistoParameters(fname, xmin, xmax);
	}
	flist.close();
	flist.open( flistname.c_str(), ios::in);
	// double dZ = (xmax-xmin)/(0.5*nbins);
    double dZ = (xmax-xmin)/(nbins);

	double k,l;
	const int MAX_DATA_POINTS(1000);
    double x, y;
	double* x_array = new double[MAX_DATA_POINTS];
	double* y_array = new double[MAX_DATA_POINTS];
	TH1D* h1 = new TH1D("h1", " ", nbins, xmin-0.5*dZ, xmax+0.5*dZ);
	
	int fcount(0);
 	int numpoints(0);
	while ( !flist.eof() )
	{
		flist >> fname;
		cout << "fname: " << fname << endl;
		if ( !flist.eof() )
		{
			ifstream fdata( fname.c_str(), ios::in);
	
	    	// Read first the header lines of the data (i.e. all lines that start with "#")
			std::string aDataLine;
   		 	bool isHeader(true);
			while (isHeader == true)
			{
				std::getline ( fdata, aDataLine);
		        isHeader = IsHeader(aDataLine);
			}
    		numpoints = 0;
		    while ( !fdata.eof() )
			{
				fdata >> x >> y;
		        if ( !fdata.eof() ) 
		        {
		            x_array[numpoints] = x;
					if (fcount == 0)
			            y_array[numpoints] = y; 
					else
			            y_array[numpoints] += y; 
					h1->Fill(x, y);
		            numpoints++;
		        }
			}
			fdata.close();
			fcount++;
		}	// if ( !flist.eof() )
	} // while ( !flist.eof() )

	// Esthetical stuff to make the plot look nice
	TStyle *style = new TStyle("gdlPRL","dry style for PRL");
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetFrameBorderMode(0);
	style->SetFrameFillStyle(0);
	style->SetCanvasColor(0);
	style->SetPadColor(0);
	style->SetOptStat(0);
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

    // This is where the data (x_array and y_array) are put into a "TGraph" (of ROOT)
	TGraph* gl4 = new TGraph(numpoints, x_array, y_array);

	gl4->SetLineColor(kBlue);
	gl4->SetLineWidth(2);
	gl4->GetXaxis()->SetTitle("Distance (mm)");

    // This is where the TGraph is drawn.
    // We don't use a histogram in this case, so we cannot adjust the bounds of the final image
    // If you want to zoom, you can:
    //  1) delete the lines of the data that you don't want to plot (from the data-file)
    //  2) change the code so that it doesn't read all lines  (but that's tricky)
    //  the easiest way: 3) do it by hand when you have ROOT open
    //
	gl4->Draw("AL");

	cout << "Do you want to Fit? (1/0)" << endl;
	bool drawFit;
	cin >> drawFit;
	if (drawFit)
	{
		gStyle->SetOptFit(1);
		gl4->Fit("gaus");
	}

	TCanvas* c2 = new TCanvas("c2","c2", 600, 600);
	c2->SetGridx();
	c2->SetGridy();
	h1->Draw("HIST");

	return 0;
}

// ============================================

bool IsHeader(const string& in_Line)
{
    const string searchStr("#");
    size_t position = in_Line.find( searchStr );
    if (position == string::npos)
        return false;
    else
        return true;
}

// ============================================

int GetHistoParameters(const string& fname, double& xmin, double& xmax)
{
	ifstream fdata( fname.c_str(), ios::in);
   	// Read first the header lines of the data (i.e. all lines that start with "#")
	std::string aDataLine;
   	bool isHeader(true);
	while (isHeader == true)
	{
		std::getline ( fdata, aDataLine);
        isHeader = IsHeader(aDataLine);
	}
   	int nbins = 0;
	double x, y;
	while ( !fdata.eof() )
	{
		fdata >> x >> y;
	    if ( !fdata.eof() ) 
	    {
			if (nbins == 0)
				xmin = x;
			xmax = x;
		    nbins++;
	    }
	}
	fdata.close();
	return nbins;
}

