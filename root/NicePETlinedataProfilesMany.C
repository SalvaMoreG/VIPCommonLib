
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

void MakeLineProfilePlot( const std::string& in_fname );
bool IsHeader(const string& aLine);
void PrintCanvas( const std::string in_fname, TCanvas* in_c );
void SetROOTStyle();

bool m_normalizepeakto1( false );

// =====================================================

int Main()
{
    std::string fname;
    cout << "give name of tsv file containing 1 or more lineprofiles" << endl;
    cin >> fname;

    cout << "dealing with file: " << fname << endl;
    MakeLineProfilePlot( fname );

    return 0;
}

// =====================================================

void
MakeLineProfilePlot( const std::string& in_fname )
{
    // Define all the style stuff (only once)
    SetROOTStyle();

    // Define the canvas
    TCanvas* c1 = new TCanvas("c1","c1", 1200, 1000);
    c1->SetGridx();
    c1->SetGridy();

    c1->SetBorderSize(10);
    c1->SetFrameBorderSize(10);
    c1->SetLineWidth(10);
    c1->SetFrameLineWidth(3);         // <---- it's this one, the others are useless (so why do they exist?)

    TLegend* legend = new TLegend(0.6,0.65,0.99,0.90);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.02);

    TCanvas* c2 = new TCanvas("c2", " ", 1000, 1000);
    c2->SetGridx();
    c2->SetGridy();

    // FILE
    ifstream fdata;
    fdata.open( in_fname.c_str(), ios::in);
    if (!fdata.is_open())
    {
        cout << "ERROR: File not open: " << in_fname << endl;
    }

    // LOOP OVER THE ENTIRE FILE TO GET ALL PROFILE NAMES
    const int maxLegend(10);
    int colors[maxLegend] = {kRed, kBlack, kBlue, kGreen+3, kMagenta+3, 
                             kOrange+7, kMagenta, kCyan, kYellow+2, kViolet+6 };
    string titles[maxLegend];
    int numTitles(0);
    string dumStr;
    string dumStr2;
    bool nextYes(false);
    while ( !fdata.eof() && numTitles < maxLegend)
    {
        fdata >> dumStr;
        if (dumStr == "Profile")
        {
            fdata >> dumStr;
            if (dumStr == "on:")
            {
                fdata >> dumStr;
				cout << "Line-profile name: " << dumStr << " Give legend tag: " << endl;
				cin >> dumStr;
                titles[numTitles] = dumStr;
                numTitles++;
            }
        }
    }

    // close and open again
    fdata.close();
    fdata.open( in_fname.c_str(), ios::in);

    // LOOP OVER THE ENTIRE FILE
    int loop(0);
    while ( !fdata.eof() )
    {
        // First read the header
        std::string aLine;
        bool isHeader(true);
        streampos oldpos;
        int numHlines( 0 );
        while (isHeader && !fdata.eof())
        {
            oldpos = fdata.tellg();  // stores the position
            fdata >> aLine;
            isHeader = IsHeader(aLine);
            if ( !isHeader )
                fdata.seekg (oldpos);
            else
            {
                fdata.ignore(1024, '\n');
                numHlines++;
            }
        }
        if (fdata.eof())
            break;

        // Read the data
        int numDlines = 0;
        double* x = new double[500];
        double* y = new double[500];
        double tmpx,tmpy;
        double maxval;
        while (!fdata.eof() && !isHeader)
        {
            oldpos = fdata.tellg();  // stores the position
            fdata >> aLine;
            isHeader = IsHeader(aLine);

            if (fdata.eof())
                break;

            if (!fdata.eof())
                fdata.seekg (oldpos);

            // Reading the actual data
            if (!isHeader && !fdata.eof())
            {
                fdata >> x[numDlines] >> y[numDlines];  // reading real values!
                if (numDlines == 0)
                {
                    maxval = y[numDlines];
                }
                else if ( y[numDlines] > maxval )
                {
                    maxval = y[numDlines];
                }
            }

            if (!fdata.eof() && !isHeader)
                numDlines++;
        }
        cout << "last data: " << x[numDlines-1] << " " << y[numDlines-1] << endl;
        cout << "number of points: " << numDlines << endl;

        double binsize = x[1] - x[0];
        double xmin = x[0] - 0.5 * binsize;
        double xmax = x[numDlines-1] + 0.5 * binsize;
        TH1D* h_tmp = new TH1D("h_tmp", " ", numDlines, xmin, xmax);
        for (int i = 0; i < numDlines; i++)
            h_tmp->Fill(x[i], y[i]);

        // Define the current graph
        if (m_normalizepeakto1 && maxval > 0.0)
        {
            for (int isigh = 0; isigh < numDlines; isigh++)
            {
                y[isigh] = y[isigh]/maxval;
            }
        }
        TGraph* gl4 = new TGraph(numDlines, x, y);

        gl4->SetTitle(" ");
        int color = loop;
        if (loop < maxLegend)
            color = colors[loop]; 
        gl4->SetLineColor(color);
        gl4->SetLineWidth(3);
        gl4->GetXaxis()->SetTitle("Distance (mm)");
        h_tmp->SetLineColor(color);

        // Draw the graphs
        if (loop == 0)
        {
            c1->cd();
            double gmax(maxval);
            gmax = gmax * 1.5;
            if (gmax > 0) 
			{
                gl4->SetMaximum(gmax);
				h_tmp->SetMaximum(gmax);
			}
            gl4->DrawClone("AL");
            c2->cd();
            h_tmp->Draw("HIST");
        }
        else
        {
            c1->cd();
            gl4->DrawClone("L SAME");
            c2->cd();
            h_tmp->Draw("HIST SAME");
        }

        string lineTitle;
        if (loop < numTitles)
            lineTitle = titles[loop];
        legend->AddEntry(gl4, lineTitle.c_str(), "L");

        loop++;
    }
    // close the file
    fdata.close();

    // Print the canvasses
    c1->cd();
    legend->DrawClone();
    PrintCanvas( in_fname, c1 );

    string histname = "hist_";
    histname += in_fname;
    cout << "histname: " << histname << endl;
    c2->cd();
    legend->DrawClone();
    PrintCanvas( histname, c2 );

    return;
}

// ==========================================================

bool IsHeader(const string& in_Line)
{
    const string searchStr("#");
    size_t position = in_Line.find( searchStr );
    if (position == string::npos)
        return false;
    else
        return true;
}

// ==========================================================

void
PrintCanvas( const std::string in_fname, TCanvas* in_c )
{
    std::string outname;
    std::string outnamepdf;
    int pos = in_fname.find(string(".tsv"));
    int len = pos;
    if (pos != std::string::npos)
    {
        outname = in_fname.substr ( 0, len );
        outnamepdf = outname;
        outname += ".gif";
        outnamepdf += ".pdf";
    }
    else
    {
        outname = in_fname;
        outnamepdf = outname;
    }

    in_c->Print( outname.c_str(), "gif");
    // in_c->Print( outnamepdf.c_str(), "pdf");
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

// ===============================================================

void SetROOTStyle()
{
    TStyle *style = new TStyle("gdlPRL","dry style for PRL");
    style->SetCanvasBorderMode(0);
    style->SetCanvasBorderSize(3);
    style->SetPadBorderMode(0);
    style->SetFrameBorderMode(0);
    style->SetFrameBorderSize(3);

    style->SetFrameFillStyle(0);
    style->SetCanvasColor(0);

    style->SetPadColor(0);
    style->SetOptStat(0);
    style->SetTextFont(132);

    style->SetPadLeftMargin(10);

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
}










