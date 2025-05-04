#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

bool IsHeader(const string& in_Line);
int GetHistoParameters(const string& fname, double& xmin, double& xmax);

enum ALGOTYPE
{
	ALG_NONE = 0,
	ALG_OSEM,
	ALG_FBP, 
	ALG_OE
};

// ============================================

int Main()
{
	ALGOTYPE algoType = ALG_NONE;
	
	double k,l;
	ifstream fdata;
	string fname("");
	cout << "give tsv filename" << endl;
	cin >> fname;
    fdata.open( fname.c_str(), ios::in);
    if (!fdata.is_open())
    {
        cout << "ERROR: File not open: " << fname << endl;
        return 1;
    }
    fdata.close();

    // Read first the header lines of the data
    double xmin, xmax;
	int nbins = GetHistoParameters(fname, xmin, xmax);

	const int MAX_DATA_POINTS(1000);
    double x, y;
	double* x_array = new double[MAX_DATA_POINTS];
	double* y_array = new double[MAX_DATA_POINTS];
    int numpoints(0);

    // double dZ = (xmax-xmin)/(0.5*nbins);
    double dZ = (xmax-xmin)/(nbins);
    TH1D* h1 = new TH1D("h1", " ", nbins, xmin-0.5*dZ, xmax+0.5*dZ);

    fdata.open( fname.c_str(), ios::in);
    std::string aDataLine;
    bool isHeader(true);
    while (isHeader == true)
    {
        std::getline ( fdata, aDataLine);
        isHeader = IsHeader(aDataLine);
		// if (isHeader) cout << "header: " << aDataLine << endl;
    }
	double sumcontents(0);
    while ( !fdata.eof() )  // while not at the end of file (eof)
	{
		fdata >> x >> y;
        if ( !fdata.eof() ) // if we are not at the end of the file, the x and y can be used
        {
            x_array[numpoints] = x;
            y_array[numpoints] = y;
            h1->Fill(x, y);
            numpoints++;
			sumcontents += y;
        }
	}
	cout << "Total histogram entries/contents: " << sumcontents << endl;
    fdata.close();

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

	if (algoType == ALG_NONE)
		gl4->SetTitle("");
	else if (algoType == ALG_OSEM)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. OSEM");
	else if (algoType == ALG_FBP)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. FBP");
	else if (algoType == ALG_OE)
		gl4->SetTitle("Derenzo, 1.02mm rods line profile. OE");

	gl4->SetLineColor(kBlue);
	gl4->SetLineWidth(2);
	gl4->GetXaxis()->SetTitle("Distance (mm)");

	if (algoType == ALG_OSEM)
		gl4->SetMaximum(0.2);	// OSEM
	else if (algoType == ALG_FBP)
		gl4->SetMaximum(200);	// FBP 
	else if (algoType == ALG_OE)
		gl4->SetMaximum(80000);	// OE 

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

    bool drawlines(false);
    // ALL THE STUFF BELOW IS ONLY IF YOU WANT TO DRAW HORIZONTAL LINES THROUGH THE HISTOGRAM, AT CERTAIN HEIGHTS
    // NORMALLY WE DON'T NEED THIS.
    // (So I commented out the question)

    // cout << "Do you want to draw horizontal lines? (1/0)" << endl;
    // cin >> drawlines;

	/*
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
        lineAve.SetX1(x_array[0]);
        lineAve.SetX2(x_array[n-1]);
        lineAve.SetY1(value);
        lineAve.SetY2(value);
        lineAve.DrawClone();

        SNR = value;

        argh = int (10 * value + 0.5);
        sigh = argh/10.0;

        TString str1("<peak>: ");
        str1 += sigh;
        lineText.SetText(x_array[0], value, str1 );
        lineText.DrawClone();

        double SNRypos = 0.5 * value;

        cout << "give <valley>" << endl;
        cin >> value;
        lineAve.SetX1(x_array[0]);
        lineAve.SetX2(x_array[n-1]);
        lineAve.SetY1(value);
        lineAve.SetY2(value);
        lineAve.DrawClone();

        SNR = SNR/value;

        argh = int (10 * value + 0.5);
        sigh = argh/10.0;
        TString str2("<valley>: ");
        str2 += sigh;
        lineText.SetText(x_array[0], value, str2 );
        lineText.DrawClone();

        SNR = (int (10 * SNR) + 0.5)/10.0;
        TString str3("SNR: ");
        str3 += SNR;
        lineText.SetText(x_array[0], SNRypos, str3);
        lineText.DrawClone();
    }
	*/

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


