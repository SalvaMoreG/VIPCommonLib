
#include <iostream>
#include <fstream>
#include <random>
    //  C++ random generator

double GetRandomUniformValue(const double& in_min, const double& in_max);

using namespace std;

const double kPI = 3.14159265358979323846;
const double PET_radius(120.0);

#define USE_ROOT 1

bool debug(false);

#ifdef USE_ROOT
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLine.h"
#include "TStyle.h"
#endif

// ==================================================

int main()
{
    // Parameters
    double source_x(0.0), source_y(0.0), source_z(0.0);
    double SIZE_SRC_XY(20.0), SIZE_SRC_Z(0.0);
    double MAX_HIT_Z(150.0);    

    ofstream outfile("COINC_PETLORS.dat_NEW");
    // ofstream outsrcfile("SOURCEPOS.dat_NEW");

    // RANDOM STUFF    
    double angleMin = 0.0; double angleMax = 2.0*kPI;
	    
    double thetaRad(0.0);
    double x1(0.0), y1(0.0), z1(0.0);
    double x2(0.0), y2(0.0), z2(0.0);    
    
    cout << "Give source size for X and Y" << endl;
    cout << "only give one number, the size is the same for X and Y" << endl;
    cin >> SIZE_SRC_XY;
        
    bool is2D(false);
    cout << "is 2D?" << endl;
    cin >> is2D;
    if (is2D)
        MAX_HIT_Z = 0.0;
    else
    {
        cout << "Give source size for Z (multiple of 10 mm)" << endl;
        cin >> SIZE_SRC_Z;
    }
    
#ifdef USE_ROOT  
    gStyle->SetOptStat(0);
    // a 2D histogram to show the position of all hits
    TH2D* h1 = new TH2D("h1", "HITS", 2*PET_radius+1, -PET_radius, PET_radius
                                    , 2*PET_radius+1, -PET_radius, PET_radius);
    
    // a 2D histogram of twice the size of the source CUBE, 
    //  to show the source position distribution, with pixel size = 2
    TH2D* h2 = new TH2D("h2", "SOURCE", SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY
                                      , SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY);
    
    // a 3D histogram of twice the size of the source CUBE, 
    //  to show the source position distribution, with pixel size (XY) = 2 and pixel size (Z) = 10
    TH3D* h3 = 0;
    if (!is2D)
    {
        double histZlength = SIZE_SRC_Z + 20; // extra 10 mm on both sides
        int nbinsZ = (int) (histZlength/10.0);
        h3 = new TH3D("h3", "SOURCE-3D", SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY
                                      , SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY 
                                      , nbinsZ     , -0.5 * histZlength,  0.5 * histZlength );
        
        cout << "SIZE_SRC_Z: " << SIZE_SRC_Z << endl;
        cout << "histZlength:" << histZlength << endl;
        cout << "nbinsZ: " << nbinsZ << endl;
        cout << " h3 Z bins: " << h3->GetNbinsZ() << endl;
        
        
    }
    const int MAX_LINES(500);
    
    TLine* lines[MAX_LINES];
    TLine* linesXZ[MAX_LINES];
    int nLinesDone(0);
#endif

    cout << "Give number of LORs" << endl;
    int n_LORs; cin >> n_LORs;

    for (int iloop = 0; iloop < n_LORs; iloop++)
    {
        // Get a random point inside a rectangle in the XY plane
        //      X = [-0.5 * size; 0.5 * size]
        //      Y = [-0.5 * size; 0.5 * size]
        //      Z = 0
        source_x = GetRandomUniformValue(0, SIZE_SRC_XY) - 0.5 * SIZE_SRC_XY;
        source_y = GetRandomUniformValue(0, SIZE_SRC_XY) - 0.5 * SIZE_SRC_XY;
        
        // Z position of source
        source_z = 0.0;
        if (!is2D)
            source_z = GetRandomUniformValue(0.0, SIZE_SRC_Z)- 0.5 * SIZE_SRC_Z;
                    
        // Get an angle for the LOR
		thetaRad = GetRandomUniformValue(angleMin, angleMax);
        if (debug) cout << "theta: " << thetaRad << endl;
        
        // We just create a distance dX and dY away from the source
        // 
        // NOTE: The hits should lie outside the size of the final image. 
        // If we have a radius R (along the x-axis), 
        //    a distance S from source to the center of the PET ring
        //    and an unknown distance D from the source to the hit, then: 
        //  S^2 + D^2 - 2 * S*D*cos(angle_SD) = R^2
        //  => D^2 - 2 * S*D*cos(angle_SD) + (S^2 - R^2) = 0. 
        // This is an quadratic equation. 
        // To avoid this, we do something simpler: 
        // D = the distance D from a position inside the source region to hit 1 or hit2
        //  Just make sure D >= R + half the size (in XY) of the final image
        //  We want to make an image of twice the size of the source. So, half the size = SIZE_SRC_XY
        // 
        
        double dX = (PET_radius + SIZE_SRC_XY) * cos(thetaRad);
        double dY = (PET_radius + SIZE_SRC_XY) * sin(thetaRad);
        
        // hit 1
        x1 = source_x + dX;
        y1 = source_y + dY;
        z1 = 0.0;
        
        // hit 2
        x2 = source_x - dX;
        y2 = source_y - dY;
        z2 = 0.0;
        
        // Z positions of the hits...
        //
        if (!is2D)
        {
            double delZ = 2.0*GetRandomUniformValue(0.0, 1.0) - 1.0;
            double maxD_srchit(0.0);
            if (source_z >= 0)
                maxD_srchit = MAX_HIT_Z - source_z;
            else if (source_z < 0)
                maxD_srchit = source_z + MAX_HIT_Z;
            delZ = delZ * maxD_srchit;  
            
            z1 = source_z + delZ;
            z2 = source_z - delZ;
        }
        
        outfile << z1 << " " << y1 << " " << x1 << " " << 511.0 << "   "
                << z2 << " " << y2 << " " << x2 << " " << 511.0 << endl;
        
#ifdef USE_ROOT                
        if (iloop < MAX_LINES)
        {
            h1->Fill(x1, y1);
            h1->Fill(x2, y2);

            TLine* aLine = new TLine(x1, y1, x2, y2);
            lines[iloop] = aLine;
            
            nLinesDone++;
        }
        h2->Fill(source_x, source_y);
        if (!is2D)
            h3->Fill(source_x, source_y, source_z);
#endif
    }
    
#ifdef USE_ROOT    
// SOURCE CUBE DENSITY PER PIXEL
    if (is2D)
    {
        ofstream binaryfile("SOURCEPOS_AMIDE.bin_NEW", ios::binary | ios::out );
        for (int biny = 1; biny <= h2->GetNbinsY(); biny++) // ROOT starts counting with 1
        {
            for (int binx = 1; binx <= h2->GetNbinsX(); binx++) // ROOT starts counting with 1
            {
                double bcx = ((TAxis*)h2->GetXaxis())->GetBinCenter(binx);
                double bcy = ((TAxis*)h2->GetYaxis())->GetBinCenter(biny);            
                float value = h2->GetBinContent (binx, biny);
                binaryfile.write((char*) &value, sizeof(value));
            }
        }
        binaryfile.close();
        ofstream headerfile("2Dheader.hv");
        headerfile << "BINARY 2D FILE: " << h2->GetNbinsX() 
            << " * " << h2->GetNbinsY() << endl;
        headerfile.close();        
    }
    else
    {
        ofstream binaryfile("SOURCEPOS_3D_AMIDE.bin_NEW", ios::binary | ios::out );
        for (int binz = 1; binz <= h3->GetNbinsZ(); binz++) // ROOT starts counting with 1
        {    
            for (int biny = 1; biny <= h3->GetNbinsY(); biny++) // ROOT starts counting with 1
            {
                for (int binx = 1; binx <= h3->GetNbinsX(); binx++) // ROOT starts counting with 1
                {
                    double bcx = ((TAxis*)h3->GetXaxis())->GetBinCenter(binx);
                    double bcy = ((TAxis*)h3->GetYaxis())->GetBinCenter(biny);
                    double bcz = ((TAxis*)h3->GetZaxis())->GetBinCenter(binz);
                    float value = h3->GetBinContent (binx, biny, binz);     
                    binaryfile.write((char*) &value, sizeof(value));
                }
            }
        }
        binaryfile.close();    
        ofstream headerfile("3Dheader.hv");
        headerfile << "BINARY 3D FILE: " << h3->GetNbinsX() 
            << " * " << h3->GetNbinsY() << " * " << h3->GetNbinsZ() << endl;
        headerfile.close();
    }

    TCanvas* m1 = new TCanvas("m1", "HITS", 700, 700);
    h1->SetMarkerStyle(20);
    h1->Draw("P");
    for (int iloop = 0; iloop < nLinesDone; iloop++)
    {
        TLine* aLine = lines[iloop];
        aLine->SetLineColor(kBlue);
        aLine->Draw("SAME");
    }
    m1->Print("hits.gif", "gif");
    
    TCanvas* m2 = new TCanvas("m2", "SOURCE", 700, 700);
    h2->Draw("COLZ");
    m2->Print("source.gif", "gif");    
#endif

	return 0;
}

// ===================

double GetRandomUniformValue(const double& in_min, const double& in_max)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(in_min, in_max);
    double randomX = dist(mt);
    if (debug)
    {
        cout << "CHECK RANDOM: " << randomX 
             << " range: " << in_min << " " << in_max << endl;
    }
    return randomX;
}



