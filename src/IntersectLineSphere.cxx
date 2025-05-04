
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>

#include "../include/CVIP3Vector.h"
#include "../include/CVIPUtils.h"

using namespace std;
using namespace VIPUtils;

void TestNEvents();
void ReadPETLors( const string& in_fname );

const double mu_value_material = 0.00968;
ofstream outfile("COINC_D_and_ATTCORR.dat_NEW");

// ==========================================

int main()
{
    string fname;
    cout << "give COINC filename" << endl;
    cout << " (if filename = 0, we are just gonna create N events): " << endl; 
    cin >> fname;
    
    if (fname.size() == 0)
        TestNEvents();
    else
        ReadPETLors( fname );
}

// =================================================================================

void TestNEvents()
{
    C3Vector X1(-120, 0.0, 0.0);
    C3Vector X2(120, 0.0, 0.0);
    
    int nEvents(0);
    cout << "How many events?" << endl;
    cin >> nEvents;
    
    const C3Vector Centre(0.0, 0.0, 0.0);
    const double radius(50.0);

	for (int ievent = 0; ievent < nEvents; ievent++)
	{
        double distance;
        bool ok = CalculateIntersectionDistance(X1, X2, Centre, radius, distance);
        if (ok)
        {
            double attcor = exp(-1.0 * mu_value_material * distance);
            outfile << X1.GetZ() << " " << X1.GetY() << " " << X1.GetX() 
                    << "   " << 511 << "     "
                    << X2.GetZ() << " " << X2.GetY() << " " << X2.GetX() 
                    << "   " << 511 << "     "
                    << distance << "  " << attcor << endl;
        }             
		X1 = X1 + C3Vector(0.0, 5.0, 0.0);
        X2 = X2 + C3Vector(0.0, 5.0, 0.0);
	}	
}

// =================================================================================

void ReadPETLors( const string& in_fname )
{
    ifstream infile(in_fname.c_str());
    
    double x1, y1, z1, E1, x2, y2, z2, E2;
    const C3Vector Centre_source(0.0, 0.0, 0.0);
    const double radius_source(50.0);
    bool debug(false);
    
    while ( !infile.eof() )
    {
        infile >> z1 >> y1 >> x1 >> E1
               >> z2 >> y2 >> x2 >> E2;
        infile.ignore(1024, '\n');
        if ( !infile.eof() )
        {
            C3Vector X1(x1, y1, z1);
            C3Vector X2(x2, y2, z2);
            
            double distance;
            bool ok = CalculateIntersectionDistance(X1, X2, Centre_source
                , radius_source, distance);
            
            if (ok)
            {
                double attcor = exp(-1.0 * mu_value_material * distance);
                if (debug)
                {
                    cout << "DEBUG: <<<<<<<" << endl;
                    cout << "X1 " << X1 << " X2: " << X2 << endl;
                    cout << "distance: " << distance 
                    << " attcor: " << attcor << endl;
                    cout << " >>>>>>>> " << endl;
                    
                    int dummy(0); cin >> dummy;
                }
                outfile << X1.GetZ() << " " << X1.GetY() << " " << X1.GetX() 
                        << "   " << 511 << "     "
                        << X2.GetZ() << " " << X2.GetY() << " " << X2.GetX() 
                        << "   " << 511 << "     "
                        << distance << "  " << attcor << endl;
            }
        }
    }
}

// =================================================================================



	
	
