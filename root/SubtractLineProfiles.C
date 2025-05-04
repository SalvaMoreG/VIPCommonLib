#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

bool IsHeader(const string& in_Line);
bool ReadFile(const string& in_fname, double out_x_array[], double out_y_array[],
	int& out_numpoints, int& out_totcounts);

const int MAX_DATA_POINTS(500);

int Main()
{
	double x_array_1[MAX_DATA_POINTS];
	double y_array_1[MAX_DATA_POINTS];
	int numpoints_1(0), totcounts_1(0);
	double x_array_2[MAX_DATA_POINTS];
	double y_array_2[MAX_DATA_POINTS];
	int numpoints_2(0), totcounts_2(0);

	string fname("");
	cout << "give filename 1" << endl;
	cin >> fname;
	bool ok = ReadFile(fname, x_array_1, y_array_1, numpoints_1, totcounts_1);
	if (!ok)
	{
		cout << "ERROR reading file: " << fname << endl;
		return 1;
	}
	cout << "fname: " << numpoints_1 << " tot contents: " << totcounts_1 << endl;
	cout << "give filename 2" << endl;
	cin >> fname;
	ok = ReadFile(fname, x_array_2, y_array_2, numpoints_2, totcounts_2);
	if (!ok)
	{
		cout << "ERROR reading file: " << fname << endl;
		return 1;
	}
	cout << "fname: " << numpoints_2 << " tot contents: " << totcounts_2 << endl;

	if (numpoints_1 != numpoints_2)
	{
		cout << "ERROR! #points in file-1 != #points in file-2 " << endl;
		cout << "Make sure the line-profiles were made from images with the same FOV!" << endl;
		cout << "CC_OSEM should have used the same osem_fov_parameters.conf" << endl;
		return 1;
	}

	// Canvas
	TCanvas* m_1 = new TCanvas("m_1", " ", 1000, 1000);
	m_1->Divide(2, 2);
    // This is where the data (x_array and y_array) are put into a "TGraph" (of ROOT)
	m_1->cd(1);
	TGraph* gl4_1 = new TGraph(numpoints_1, x_array_1, y_array_1);
	gl4_1->SetLineColor(kBlue);
	gl4_1->SetLineWidth(2);
	gl4_1->GetXaxis()->SetTitle("Distance (mm)");
	gl4_1->Draw("AL");

	m_1->cd(2);
	TGraph* gl4_2 = new TGraph(numpoints_2, x_array_2, y_array_2);
	gl4_2->SetLineColor(kRed);
	gl4_2->SetLineWidth(2);
	gl4_2->GetXaxis()->SetTitle("Distance (mm)");
	gl4_2->Draw("AL");

	// Normalize
	double factor = (double) (totcounts_1)/(double) (totcounts_2);
	double y_array_norm_2[MAX_DATA_POINTS];
	double y_subtract_1_min_2[MAX_DATA_POINTS];
	for (int i = 0; i < numpoints_2; i++)
	{
		y_array_norm_2[i] = factor*y_array_2[i];
		y_subtract_1_min_2[i] = y_array_1[i] - y_array_2[i];
	}
	
	m_1->cd(3);
	gl4_1->Draw("AL");
	TGraph* gl4_3 = new TGraph(numpoints_2, x_array_2, y_array_norm_2);
	gl4_3->SetLineColor(kRed);
	gl4_3->SetLineWidth(2);
	gl4_3->Draw("L SAME");

	m_1->cd(4);
	TGraph* gl4_4 = new TGraph(numpoints_2, x_array_2, y_subtract_1_min_2);
	gl4_4->SetLineColor(kGreen+4);
	gl4_4->SetLineWidth(2);
	gl4_4->Draw("AL");

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

bool ReadFile(const string& in_fname, double out_x_array[], double out_y_array[],
	int& out_numpoints, int& out_totcounts)
{
	ifstream fdata;
	fdata.open( in_fname.c_str(), ios::in);
	if (!fdata.is_open()) 
	{
		cout << "ERROR: File not open: " << in_fname << endl;
		return false;
	}

	std::string aDataLine;
    bool isHeader(true);
	while (isHeader == true)
	{
		std::getline ( fdata, aDataLine);
        isHeader = IsHeader(aDataLine);
	}

    double x, y;
    out_numpoints = 0;
	out_totcounts = 0;
    while ( !fdata.eof() )  // while not at the end of file (eof)
	{
		fdata >> x >> y;
        if ( !fdata.eof() ) // if we are not at the end of the file, the x and y can be used
        {
            out_x_array[out_numpoints] = x;
            out_y_array[out_numpoints] = y;
            out_numpoints++;
			out_totcounts += y;
        }
	}
	return true;
}









