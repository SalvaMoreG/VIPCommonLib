
#include "../include/CVIPToyPET_aux.h"

#include "../include/CVIP3Vector.h"
#include "../include/VIPconstants.h"

#include "../include/CVIPRandom.h"

#include <iostream>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLine.h"
#include "TStyle.h"

using namespace std;

void CalculateHit1ToSource( const C3Vector& in_source, const C3Vector& in_hit1, C3Vector& out_hit1Tosrc );

void CalculateSourceToOrigin( const C3Vector& in_source, C3Vector& out_source2O );

double CalculateCosAngleBD( const C3Vector& in_directionD, const C3Vector& in_directionB );

//  int CalculateLengthBOfTriangle(const double& in_a, const double& in_b, const double& in_cosAngle,
//                             double& out_solution1, double& out_solution2 );

void CalculateHit2( const C3Vector& in_source, const C3Vector& in_direction
    , const double& in_lengthSrc2Hit2, C3Vector& out_hit2 );

// ****************************************************************
// Explanation of this code: See: CalculateLengthBOfTriangle(...)
// ****************************************************************

const double PET_radius(120.0);
const int MAX_LINES(10000);

int main()
{
    gStyle->SetOptStat(0);

    // Parameters
    double source_x(0.0), source_y(0.0), source_z(0.0);
    
    double SRC_CTR_X(0.0), SRC_CTR_Y(0.0);
    double SIZE_SRC_XY(80.0), SIZE_SRC_Z(0.0);
    double MAX_HIT_Z(150.0);
    
    double PET_voxelSizeZ(0.0);    

    // cout << "Fixed source (0) or Random source (!=0)?" << endl;
    bool doRandom(true);
    bool is2D(false);
    // cin >> doRandom;
    if (doRandom == 0)
    {
        cout << "Give source x and y < PET_radius: " << PET_radius << endl;
        cin >> source_x >> source_y;
    }
    else
    {
        cout << "Give source size for X and Y" << endl;
        cout << "only give one number, the size is the same for X and Y" << endl;
        cin >> SIZE_SRC_XY;
        
        cout << "is 2D?" << endl;
        cin >> is2D;
        
        if (is2D) 
        {
            MAX_HIT_Z = 0.0;
        }
        else
        {
            cout << "Give source size for Z" << endl;
            cin >> SIZE_SRC_Z;
            
            cout << "Give PET Voxel Size in Z >= 0.0 && <= 10 mm (0 = use real impact point)" 
                 << endl;
            cin >> PET_voxelSizeZ; 
            if (PET_voxelSizeZ < 0 || PET_voxelSizeZ > 10)
                PET_voxelSizeZ = 0.0;
            cout << "PET_voxelSizeZ: " << PET_voxelSizeZ << endl;
            
        }
    }
    C3Vector sourcePos( source_x, source_y, 0.0);
    
    ofstream outfile("COINC_PETLORS.dat_NEW");

    // RANDOM STUFF
    CVIPRandomUniform uniformSrcDist(0.0, SIZE_SRC_XY);
    CVIPRandomUniform uniformSrcDistZ(0.0, SIZE_SRC_Z);
    
	CVIPRandomUniform uniformHitDist(0.0, 2.0*kPI);
    CVIPRandomUniform uniformHitDistZ(0.0, 1.0);
    
    double theta(0.0), x(0.0), y(0.0);
    C3Vector hit1, hit2;
    C3Vector hit1ToSrc, srcToO;
    double lengthB, cosAngle, solution1, solution2;

    TH2D* h1 = new TH2D("h1", " ", 2*PET_radius+1, -PET_radius, PET_radius
                                 , 2*PET_radius+1, -PET_radius, PET_radius);
    /*
    TH2D* hxz = new TH2D("hxz", " ", 2*PET_radius+1, -PET_radius, PET_radius
                                 , 2*PET_radius+1, -PET_radius, PET_radius);
    */
    
    // a 2D histogram of twice the size of the source CUBE, to show the source position distribution
    // , with pixel size = 2
    TH2D* h2 = new TH2D("h2", "SOURCE", SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY
                                      , SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY);

    // a 3D histogram of twice the size of the source CUBE, to show the source position distribution
    //  , with pixel size = 2 and 10 in Z
    double histZlength = SIZE_SRC_Z + 20; // extra 10 mm on both sides
    int nbinsZ = (int) (histZlength/10.0);
    TH3D* h3 = new TH3D("h3", "SOURCE-3D", SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY
                                      , SIZE_SRC_XY, -SIZE_SRC_XY, SIZE_SRC_XY 
                                      , nbinsZ     , -0.5 * histZlength,  0.5 * histZlength );

    // h1->Fill(sourcePos.GetX(), sourcePos.GetY());
    bool debug(false);

    TLine* lines[MAX_LINES];
    //  TLine* linesXZ[MAX_LINES];
    int nLinesDone(0);
        
    cout << "Give number of LORs" << endl;
    int n_LORs; cin >> n_LORs;

    for (int iloop = 0; iloop < n_LORs; iloop++)
    {
        if (doRandom)
        {
            source_x = uniformSrcDist.GetNewValue();
            source_x += SRC_CTR_X - 0.5 * SIZE_SRC_XY;
            source_y = uniformSrcDist.GetNewValue();
            source_y += SRC_CTR_X - 0.5 * SIZE_SRC_XY;

            // CUBE
            source_z = uniformSrcDistZ.GetNewValue() - 0.5 * SIZE_SRC_Z;
            
            // SQUARE
            // source_z = source_z_NULL;
            
            // if (iloop < MAX_LINES)
            //      h3->Fill(source_x, source_y, source_z);
            
            sourcePos.Set( source_x, source_y, 0.0);
            //
            // Only later will the source Z position be adapted
            //
        }
        if (debug) cout << "SOURCE: " << sourcePos << endl;
        
		theta = uniformHitDist.GetNewValue();
        if (debug) cout << "theta: " << theta << endl;
        x = PET_radius * cos(theta);
        y = PET_radius * sin(theta);
        hit1.Set(x, y, 0.0);
        if (debug) cout << "hit1: " << hit1 << endl;

        // Get Hit 2
        CalculateHit1ToSource( sourcePos, hit1, hit1ToSrc );
        if (debug) cout << "hit1ToSrc: " << hit1ToSrc << std::endl;
        CalculateSourceToOrigin( sourcePos, srcToO );
        if (debug) cout << "o2rsc: " << srcToO << " with length: " << srcToO.GetLength() << endl;
        double cosAngle = CalculateCosAngleBD( hit1ToSrc, srcToO );
        if (debug) cout << "cosAngle B-D: " << cosAngle << " so angle: " 
                        << (180/kPI) * acos(cosAngle) << endl;
        int numsolutions = CalculateLengthBOfTriangle(PET_radius, srcToO.GetLength(), cosAngle, solution1, solution2 );
        if (debug)
        {
            if (numsolutions > 0)
                cout << " solution(1): " << solution1;
            if (numsolutions > 1)
                cout << " solution(2): " << solution2;
            cout << endl;
        }

        double distanceSrc2Hit2 = 0.0;
        if (numsolutions == 1)
            distanceSrc2Hit2 = solution1;
        else if (numsolutions == 2)
        {
            if (solution1 > 0 && solution2 < 0)
                distanceSrc2Hit2 = solution1;
            else if (solution1 < 0 && solution2 > 0)
                distanceSrc2Hit2 = solution2;
            else if (solution1 > 0 && solution2 > 0)
                distanceSrc2Hit2 = (solution1 > solution2) ? solution1 : solution2;
        }

        if (distanceSrc2Hit2 > 0.0)
        {
            CalculateHit2( sourcePos, hit1ToSrc, distanceSrc2Hit2, hit2 );
            if (debug) cout << "HIT 2: " << hit2 << endl;
            
            //
            // Only now, set hit1.Z and hit2.Z
            //
            
            double delZ = 2.0*uniformHitDistZ.GetNewValue() - 1.0;
            double maxD_srchit(0.0);
            if (source_z >= 0)
                maxD_srchit = MAX_HIT_Z - source_z;
            else if (source_z < 0)
                maxD_srchit = source_z + MAX_HIT_Z;
            delZ = delZ * maxD_srchit;
            
            if (PET_voxelSizeZ > 0.1)
            {
                double Zpos1 = source_z + delZ;
                double Zpos2 = source_z - delZ;
/*                
if (source_z > 35 && Zpos1 > 145)
{
    cout << "source_z: " << source_z << endl;
    cout << "maxD: " << maxD_srchit << " delZ: " << delZ << endl;
    cout << "Zpos 1: " << Zpos1 << " 2: " << Zpos2 << endl;
}
*/
                
                int nZbin = int (Zpos1/PET_voxelSizeZ);
                Zpos1 = nZbin * PET_voxelSizeZ;
                nZbin = int (Zpos2/PET_voxelSizeZ);
                Zpos2 = nZbin * PET_voxelSizeZ;
                hit1.SetZ(Zpos1);
                hit2.SetZ(Zpos2);
            }
            else
            {
                hit1.SetZ(source_z + delZ);
                hit2.SetZ(source_z - delZ);
            }
            
            if (iloop < MAX_LINES)
            {
                h1->Fill(hit1.GetX(), hit1.GetY());
                h1->Fill(hit2.GetX(), hit2.GetY());

                TLine* aLine = new TLine(hit1.GetX(), hit1.GetY(), hit2.GetX(), hit2.GetY());
                lines[iloop] = aLine;
                
                // hxz->Fill(source_x, source_z);
                // TLine* aLineXZ = new TLine(hit1.GetX(), hit1.GetZ(), hit2.GetX(), hit2.GetZ());
                // linesXZ[iloop] = aLineXZ;
                
                nLinesDone++;
            }
            h2->Fill(source_x, source_y);
            if (!is2D)
                h3->Fill(source_x, source_y, source_z);
            
            outfile << hit1.GetZ() << " " << hit1.GetY() << " " << hit1.GetX() 
                    << " " << 511.0 << "   "
                    << hit2.GetZ() << " " << hit2.GetY() << " " << hit2.GetX() 
                    << " " << 511.0 << endl;
        }
        if (debug) cout << iloop << endl;
    }

    cout << "DONE. Now some final ROOT stuff. #lines done: " << nLinesDone << endl;
    
    TCanvas* m1 = new TCanvas("m1", " ", 700, 700);
    h1->SetMarkerStyle(20);
    h1->Draw("P");
    for (int iloop = 0; iloop < nLinesDone; iloop++)
    {
        TLine* aLine = lines[iloop];
        aLine->SetLineColor(kBlue);
        aLine->Draw("SAME");
    }
    m1->Print("hits.gif", "gif");
    
    
    // SOURCE stuff
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
    ofstream header2Dfile("2Dheader.hv");
    header2Dfile << "BINARY 2D FILE: " << h2->GetNbinsX() << " * " << h2->GetNbinsY() << endl;
    header2Dfile.close();    
    
    TCanvas* m2 = new TCanvas("m2", "SOURCE", 700, 700);
    h2->Draw("COLZ");
    m2->Print("source.gif", "gif");
    
    // SOURCE 3D stuff
    if (!is2D)
    {
        binaryfile.open("SOURCEPOS_3D_AMIDE.bin_NEW", ios::binary | ios::out );
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
    
        TCanvas* m3 = new TCanvas("m3", "SOURCE 3D", 700, 700);
        h3->Draw("");
        m3->Print("source3D.gif", "gif");    
    }

	return 0;
}

void
CalculateHit1ToSource( const C3Vector& in_source, const C3Vector& in_hit1, C3Vector& out_hit1Tosrc )
{
    out_hit1Tosrc = in_source - in_hit1;
}

void
CalculateSourceToOrigin( const C3Vector& in_source, C3Vector& out_sourceToO )
{
    out_sourceToO = in_source * -1.0 ;   // i.e. origin - source
}

double
CalculateCosAngleBD( const C3Vector& in_direction1, const C3Vector& in_direction2 )
{
    double cosAngle = in_direction1.GetScalarProductCosAngle( in_direction2 );
    return cosAngle;
}

void CalculateHit2( const C3Vector& in_source, const C3Vector& in_direction
    , const double& in_lengthSrc2Hit2, C3Vector& out_hit2 )
{
    out_hit2 = in_source + in_direction * (in_lengthSrc2Hit2/in_direction.GetLength());
}










