#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

bool IsHeader(const string& in_Line);

// ============================================

int Main()
{
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

    // Get rid of the header
    std::string aDataLine;
    bool isHeader(true);
    while (isHeader == true)
    {
        std::getline ( fdata, aDataLine);
        isHeader = IsHeader(aDataLine);
		// if (isHeader) cout << "header: " << aDataLine << endl;
    }
    
    // Read in all values;
    // 
	const int MAX_DATA_POINTS(1000);
	double* x_array = new double[MAX_DATA_POINTS];
	double* y_array = new double[MAX_DATA_POINTS];    
    double x, y;
    int numpoints(0);
    double maxValue(0);
    while ( !fdata.eof() )  // while not at the end of file (eof)
	{
		fdata >> x >> y;
        if ( !fdata.eof() ) // if we are not at the end of the file, the x and y can be used
        {
            x_array[numpoints] = x;
            y_array[numpoints] = y;
            if (y > maxValue) maxValue = y;
            numpoints++;
        }
	}
	cout << "Total bins: " << numpoints << endl;
    fdata.close();

    // Find the peaks
    // 
    std::vector<int> found_peaks;
    
    double peakThreshold = 0.20 * maxValue;

    for (int ibin = 2; ibin < numpoints-2; ibin++) // skip the edges
    {
        double value = y_array[ibin];
        if (   value > peakThreshold 
            && (value == y_array[ibin-1]) )
        {
            // value += 1; // artificially make sure that this pixel has preference over the neighbours
            int dum = 0;
        }
        else if (   value > peakThreshold
            && value > y_array[ibin-1] && value >= y_array[ibin+1]
            && value > y_array[ibin-2] && value > y_array[ibin+2] 
            && value > y_array[ibin-3] && value > y_array[ibin+3] 
           )
        {
            found_peaks.push_back(ibin);
        }
    } 
    cout << "Found " << found_peaks.size() << " peaks" << endl;
    
    // Find peaks to valley
    //
    double average_p2vs(0.0);
    std::vector<double> found_p2vs;
    for (int ip = 1; ip < found_peaks.size(); ip++)
    {
        int prev_bin = found_peaks[ip-1];
        int bin = found_peaks[ip];

        int mid_bin = 0.5 * (double) (prev_bin + bin);
            // prev_bin + (bin - prev_bin)/2 = (prev_bin + bin)/2
        double valley = y_array[mid_bin];
        if (valley < 1) valley = 1;
        double prev_value = y_array[prev_bin];
        double value = y_array[bin];
        
        double p2v = ((value + prev_value)/2.0)/valley;
        found_p2vs.push_back(p2v);
        average_p2vs += p2v;  
        
        cout << "prev peak (@ " << x_array[prev_bin] << ") = " << prev_value
             << "; this peak (@ " << x_array[bin] << ") = " << value
             << "; valley (@ " << x_array[mid_bin] << ") = " << valley
             << "; peak2valley: " << p2v << endl;
    }
    
    // average p2valley and sigma
    average_p2vs = average_p2vs/found_p2vs.size();

    double sigma_p2vs = 0.0;
    for (int ip2v = 0; ip2v < found_p2vs.size(); ip2v++)
    {
        double delta = found_p2vs[ip2v] - average_p2vs;
        sigma_p2vs += delta * delta;
    }
    sigma_p2vs = sqrt(sigma_p2vs / found_p2vs.size());  
    
    cout << "Found " << found_peaks.size() << " peaks"
         << " and " << found_p2vs.size() << " peak-to-valleys"
         << endl;
    cout << " <p2v>: " << average_p2vs << " +- " << sigma_p2vs << endl;    
    
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

