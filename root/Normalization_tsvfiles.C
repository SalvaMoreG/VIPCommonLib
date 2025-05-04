#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

bool IsHeader(const string& in_Line); 
bool ReadFile(const string& in_fname, double out_x_array[], double out_y_array[],
	int& out_numpoints, double& out_totcounts);

const int MAX_DATA_POINTS(1000);

/*
// Give N tsv files. The contents of the tsv files are normalized to the sum of the contents of the 1st tsv file.
*/

int Main()
{
	double x_array[MAX_DATA_POINTS];  //////////////array 
	double y_array[MAX_DATA_POINTS];

	int numpoints(0);
	double sum_1(0), sum_2(0), y_norm(0);

	ofstream myfile( "normalized_lineprofiles.tsv_NEW" );

	string fname("");  /////string data type variable

	cout << "Give filename (the rest will be normalized to this one)" << endl;
	cin >> fname;
	if ( !ReadFile(fname, x_array, y_array, numpoints, sum_1) )
	{
		cout << "ERROR reading file: " << fname << endl;
		return 1;
	} else {
		cout << "SUM: " << sum_1 << endl;
		myfile << "# " << endl;
		myfile << "# Profile on: " << fname << endl;
		myfile << "# " << endl;
		for(int i=0; i < numpoints; i++) 
			myfile << x_array[i] << "\t" << y_array[i] << "\n";
	}

	cout << "Give filename (contents will be normalized)" << endl;
	cin >> fname;
	while (fname.size() > 2)
	{
		if ( !ReadFile(fname, x_array, y_array, numpoints, sum_2) )
		{
			cout << "ERROR reading file: " << fname << endl;
			return 1;
		} else {
	   		cout << "SUM: " << sum_2 << endl;
		}

		double sum_ratio = (sum_1 / sum_2);
		cout << "RATIO: " << sum_ratio << endl;

		myfile << "# " << endl;
		myfile << "# Profile on: " << fname << endl;
		myfile << "# " << endl;
		for (int i=0; i < numpoints; i++) 
		{
			y_norm = (double)  (y_array[i]  * sum_ratio);
			myfile << x_array[i] << "\t" << y_norm << endl;
		}

		cout << "Give filename (contents will be normalized)" << endl;
		cin >> fname;
	}

	myfile.close();
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
	int& out_numpoints, double& out_totcounts)
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
		std::getline (fdata, aDataLine);
        isHeader = IsHeader(aDataLine);
	}

	// Add code to process the first non-header line here
	// else it will get discarded while reading

    double x, y;
    out_numpoints = 0;
	out_totcounts = 0;
    while ( !fdata.eof() )  // while not at the end of file (eof)
	{
		fdata >> x >> y;  // method for reading space separated values
        if ( !fdata.eof() ) // if we are not at the end of the file, the x and y can be used
        {
            out_x_array[out_numpoints] = x;  // Storing the values on x and y in the array with increased index at every run
            out_y_array[out_numpoints] = y;
            out_numpoints++;
			out_totcounts += y;
        }
	}
	return true;
}

