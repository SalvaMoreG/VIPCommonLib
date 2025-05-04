
#include "CVIPToyPET_aux.h"

#include "CVIP3Vector.h"
#include "CVIP3Matrix.h"
#include "VIPconstants.h"
#include "CVIPUtils.h"

#include "CVIPRandom.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TText.h"
#include "TEllipse.h"

// Main codes

int DoToyPET_3Hits_TOF();

// <<<
// Methods for just a PET ring with a random source position (see CVIPToyPET.cxx)
//  TODO: avoid code-duplication and use the functions from there
//
void CalculateHit1ToSource( const C3Vector& in_source, const C3Vector& in_hit1, C3Vector& out_hit1Tosrc );

void CalculateSourceToOrigin( const C3Vector& in_source, C3Vector& out_source2O );

double CalculateCosAngleBD( const C3Vector& in_directionD, const C3Vector& in_directionB );

// int
// CalculateLengthBOfTriangle(const double& in_A, const double& in_C, const double& in_cosAngle_A,
//                            double& out_solution1, double& out_solution2 );

void CalculateHit2( const C3Vector& in_source, const C3Vector& in_gradient, const double& in_distanceSrc2Hit2,
                    C3Vector& out_hit2 );
// >>>
// ----
// <<<
// Methods for doing the "three-hit-TOF" trick:
//

void DoThe3HitTOFStuff(int iloop, const C3Vector& in_hit1, const C3Vector& in_hit2, const C3Vector& in_hit3
    , TLine** in_lines, TLine** in_lines2nd, const C3Vector& in_sourcePos
    , const C3Vector& in_scatterPosition1, const C3Vector& in_scatterPosition2, const C3Vector& in_scatterPosition3
    , const unsigned long long int& in_eventtime_ps, CVIPRandomGauss& in_gaussDist_time
    , const C3Vector& in_orig_hit1, const C3Vector& in_orig_hit2, const C3Vector& in_orig_hit3 );

double GetSingleTimeResolutionGauss();

C3Vector CalculateSourcePosition( const C3Vector& in_hit1, const double& in_time_1,
                                  const double& in_time_2, const C3Vector& in_lor_from2_to_1 );

double GetCosAlphaThreeTriangleSides( const double& in_side1_AF, const double& in_side2_AS, const double& in_side3_FS );

double GetDeltaTimeDifference(const C3Vector& approx_source, const C3Vector& in_hit1, const double& in_time1
    , const C3Vector& in_hit3, const double& in_time3);

double UpdateApproximateSource(C3Vector& approx_source, const C3Vector& in_hit1, const double& in_time1
    , const C3Vector& in_hit3, const double& in_time3, const C3Vector& in_lor );

void GetMonteCarloKnowledgeSourcePositionFromTravellingTimes( C3Vector& out_sourceTravelling, const double& in_time_event_ps,
    const C3Vector& in_hit_LOR, const double& in_time_LORhit, const C3Vector& in_hit_3, const double& in_time_3 );

struct ToPutInPlot
{
    double time_A_clean;
    double time_B_clean;
    double time_F_clean;
    double time_A_smeared;
    double time_B_smeared;
    double time_F_smeared;
    double dT_SA_clean;
    double dT_SB_clean;
    double dT_SF_clean;
    double dT_FA_clean;
    double dT_FB_clean;
    double dT_SA_smeared;
    double dT_SB_smeared;
    double dT_SF_smeared;
    double dT_FA;
    double dT_FB;
    double length_AB;
    double so_time_AB;
    double calculated_time_SF;
    double so_length_SF;
    C3Vector Spos_AB;
    C3Vector Spos_FA;
    C3Vector Spos_FB;
};

bool MokhtarTravelAlongLOR_II(const C3Vector& in_lor_B2A, int in_direction, const C3Vector& in_lorHit
    , const C3Vector& in_hit3, const double& in_delta_T_FX_measured
    , C3Vector& io_sourceS_F );

bool MokhtarMethod_I(const double& in_timeA, const double& in_timeB, const double& in_time3
    , const C3Vector& in_hitA, const C3Vector& in_hit3, const C3Vector& in_lor, const C3Vector& in_source_AB_fromLOR
    , ToPutInPlot& io_toPutInPlot, C3Vector& io_sourceS_F);

void GetIteratedSourceApproximation( C3Vector& out_source_XF_iteration, C3Vector& in_other_hit_fromLOR
    , const C3Vector& in_source_XF_travelling, const C3Vector& in_realSourcePos
    , const C3Vector& in_hit_LOR, const double& in_time_LORhit, const C3Vector& in_hit_3, const double& in_time_3
    , TH1D* h_S_XF_srcDeviation, TH1D* h_XsF_timeDeviation, TH1D* hErrorS_AX );

// >>>

void TestFormulas();

using namespace std;

// ****************************************************************
// Explanation of this code: See: CalculateLengthBOfTriangle(...)
// ****************************************************************

// -----------------------------------------

// CONSTANTS <<<<
//
const double k_cspeed(3E-1); // mm per picoseconds

//  #define RANDOM_SOURCE 1

const double time12_uncertainty_ps_FWHM(200);        // 200 picoseconds
const double time12_uncertainty_ps_sigma( time12_uncertainty_ps_FWHM/2.355 );

const double PET_radius = 200.0; // 20.0 mm

const double E1(511.0), E2(511.0), E3(1274.5);

const int n_LORs(10000); // (10000);

const int k_plotGoodEvents(10);

// VERY UGLY global variables...
// <<<<<<<<<<<<
int ngoodevents(0);
int ngoodAorBevents(0);
int nogoodAnogoodB(0);

TH2D* h1(0);
TH1D* htsigma(0);
TH1D* hdtsigma(0);

TH1D* hErrorS_AB(0);
TH1D* hErrorS_FA(0);
TH1D* hErrorS_FB(0);
TH1D* hErrorS_averaged(0);

TH1D* hErrorS_AB_FA(0);
TH1D* hErrorS_AB_FB(0);
TH1D* hErrorS_FA_FB(0);

std::ofstream outputfile("ToyPET.out_NEW");
// >>>>>>>>>>>>

// Drawings to see what's going on
TEllipse* petRing(NULL);
TH2D* h_petRing(NULL);


// START OF MAIN
int main()
{
// double angleComptonDegrees;
// bool ok = VIPUtils::GetComptonAngleDegrees(21.0, 1274.0, angleComptonDegrees );
// if (ok)
//     cout << "angleComptonDegrees @ E1 = 14: " << angleComptonDegrees << endl;
// return 0;

    // Drawings to see what's going on
    petRing = new TEllipse(0.0, 0.0, PET_radius, PET_radius);
    h_petRing = new TH2D("h_petRing", " ", 1000, -1.1*PET_radius, 1.1*PET_radius
                                         , 1000, -1.1*PET_radius, 1.1*PET_radius);

    int retval = DoToyPET_3Hits_TOF();
    return retval;
}

// ===================================================

int DoToyPET_3Hits_TOF()
{
//      TestFormulas();
//      return 0;

// -------------------

    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    // Parameters
    bool debug(false);

    // SINGLE SOURCE
    //
    double source_Z(0.0);
    double source_x(0.0), source_y(0.0);
#ifndef RANDOM_SOURCE
    cout << "Give source x and y < 20" << endl;
    cin >> source_x >> source_y;
    C3Vector sourcePos( source_x, source_y, source_Z);
#else
    // OR: RANDOM SOURCE IN AN AREA
    cout << "Give source x and y boundaries < 10" << endl;
    cin >> source_x >> source_y;
    C3Vector sourcePos;
#endif
    cout << endl;

    // RANDOM STUFF <<<<
    //
	CVIPRandomUniform uniformDist(0.0, 2.0*kPI);
    CVIPRandomUniform uniformDist_srcX(-source_x, source_x);
    CVIPRandomUniform uniformDist_srcY(-source_y, source_y);

    // Time resolution per single time, instead of time-difference
    // t_a - t_b +- Gauss_FWHM( time12_uncertainty_ps)

    double single_t_uncertainty_sigma_ps = GetSingleTimeResolutionGauss();
    int seed(1);
    CVIPRandomGauss gaussDist_time(0.0, single_t_uncertainty_sigma_ps, seed);
    // >>

    double theta(0.0), x(0.0), y(0.0);
    C3Vector hit1;
    C3Vector hit2;
    C3Vector hit3;
    C3Vector hit1ToSrc, srcToO;
    double lengthB, cosAngle, solution1, solution2;

    // HISTOGRAMS
    h1 = new TH2D("h1", "HITS", 251, -PET_radius-10, PET_radius+10, 251, -PET_radius-10, PET_radius+10);

    htsigma = new TH1D("htsigma", "check T* - T", 100, -400, 400);
    hdtsigma = new TH1D("hdtsigma", "check (T_1 - T_2)* - (T_1 - T_2)", 100, -400, 400);

    hErrorS_AB = new TH1D("hErrorS_AB", "Distance S_REAL - S_AB(fromLOR)' (mm)", 200, -100., 100.);
    hErrorS_FA = new TH1D("hErrorS_FA", "Distance S_REAL - S_FA(Mokhtar walk)' (mm)", 200, -100., 100.);
    hErrorS_FB = new TH1D("hErrorS_FB", "Distance S_REAL - S_FB(Mokhtar walk)' (mm)", 200, -100., 100.);
    hErrorS_averaged = new TH1D("hErrorS_averaged", "Distance S_REAL - S_averaged' (mm)", 200, -100., 100.);

    hErrorS_AB_FA = new TH1D("hErrorS_AB_FA", "Distance S_REAL - (S_AB(fromLOR) && S_FA)' (mm)", 200, -100., 100.);
    hErrorS_AB_FB = new TH1D("hErrorS_AB_FB", "Distance S_REAL - (S_AB(fromLOR) && S_FB)' (mm)", 200, -100., 100.);
    hErrorS_FA_FB = new TH1D("hErrorS_FA_FB", "Distance S_REAL - (S_FA && S_FB)' (mm)", 200, -100., 100.);

#ifndef RANDOM_SOURCE
    h1->Fill(sourcePos.GetX(), sourcePos.GetY());
    if (debug) cout << "SOURCE: " << sourcePos << endl;
#endif

    bool doScatter(false), doSmearHitPos(false);
    //
    cout << "Scatter track of hit2? (1/0)" << endl;
    cin >> doScatter;
    cout << "Smear all 3 hit positions? (1/0)" << endl;
    cin >> doSmearHitPos;

    TLine* lines[n_LORs] = {0};
    TLine* lines2nd[n_LORs] = {0};

    unsigned long long int time_ps(0);
    const int delta_time_ps(1000);   // 1 ns

    // YEAH; EVENT LOOP (LOR LOOP)
    //
    for (int iloop = 0; iloop < n_LORs; iloop++)
    {
#ifdef RANDOM_SOURCE
        sourcePos.Set( uniformDist_srcX.GetNewValue(), uniformDist_srcY.GetNewValue(), source_Z);
        h1->Fill(sourcePos.GetX(), sourcePos.GetY());
        if (ngoodevents < k_plotGoodEvents)
            cout << "Source position: " << sourcePos << endl;
#endif

        // Get Hit 1
		theta = uniformDist.GetNewValue();

        if (debug) cout << "theta: " << theta << endl;
        x = PET_radius * cos(theta);
        y = PET_radius * sin(theta);
        // y = sourcePos.GetY();        <-------- DEBUGGING
        // x = sqrt(PET_radius*PET_radius - y*y);
        hit1.Set(x, y, source_Z);
        if (debug) cout << "hit1: " << hit1 << endl;

        // Get Hit 2
        CalculateHit1ToSource( sourcePos, hit1, hit1ToSrc );
        if (debug) cout << "hit1ToSrc: " << hit1ToSrc << std::endl;
        CalculateSourceToOrigin( sourcePos, srcToO );
        if (debug) cout << "o2rsc: " << srcToO << " with length: " << srcToO.GetLength() << endl;
        double cosAngle = CalculateCosAngleBD( hit1ToSrc, srcToO );
        if (debug) cout << "cosAngle B-D: " << cosAngle << " so angle: " << (180/kPI) * acos(cosAngle) << endl;
        //  int numsolutions = CalculateDistanceSrcToHit2(PET_radius, srcToO.GetLength(), cosAngle, solution1, solution2 );
        int numsolutions = CalculateLengthBOfTriangle(PET_radius, srcToO.GetLength(), cosAngle, solution1, solution2 );
            /*
            In ANY triangle with sides with lengths R, AS, SO,
            and an angle alpha_opposite_R between AS(unknown, but along known lor) and SO
            the following is true:
                AS^2 - 2 * SO * cos(alpha_opposite_R)* AS + (SO^2 - R^2) = 0
            Quadratic equation:
                unknown AS -> x
                    a = 1
                    b = - 2 * SO * cos(alpha_opposite_R)
                    c = (SO^2 - R^2)
            */

        if (debug)
        {
            if (numsolutions > 0)
                cout << " solution(1): " << solution1;
            if (numsolutions > 1)
                cout << " solution(2): " << solution2;
            cout << endl;
        }

        // Get Hit 3 (e.g. 1274 keV) <<<<
        theta = uniformDist.GetNewValue();
        if (debug) cout << "theta: " << theta << endl;
        x = PET_radius * cos(theta);
        y = PET_radius * sin(theta);
        // x = 20.0; y = 10.0;  // DEBUGGING!!!!!!
        hit3.Set(x, y, source_Z);
        // >>>>

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

        // If we could create a second hit on the ring, on the LOR that traverses the source
        if (distanceSrc2Hit2 > 0.0)
        {
            CalculateHit2( sourcePos, hit1ToSrc, distanceSrc2Hit2, hit2 );
            //  cout << "HIT 2: " << hit2 << endl;
            //  cout << "HIT 3: " << hit3 << endl;

            C3Vector orig_hit1(hit1);
            C3Vector orig_hit2(hit2);
            C3Vector orig_hit3(hit3);

            // If scatter track 2:
            C3Vector scatterPosition1(sourcePos), scatterPosition2(sourcePos), scatterPosition3(sourcePos);
            double radiusSphere(2.0);   // 2 mm
            double energy_hit1(511.0), energy_hit2(511.0), energy_hit3(1274.0);
            bool okscatter(true);

            if (doScatter)
            {
                okscatter = CalculateScatteredTrack( PET_radius, radiusSphere
                        , sourcePos, hit1, scatterPosition1, energy_hit1 );
                okscatter = okscatter && ( fabs(energy_hit1 - mass_electron_keV) < 14.0 );
                if (okscatter)
                {
                    okscatter = CalculateScatteredTrack( PET_radius, radiusSphere
                            , sourcePos, hit2, scatterPosition2, energy_hit2 );
                    okscatter = okscatter && ( fabs(energy_hit2 - mass_electron_keV) < 14.0 );
                }
                if (okscatter)
                {
                    okscatter = CalculateScatteredTrack( PET_radius, radiusSphere
                            , sourcePos, hit3, scatterPosition3, energy_hit3 );
                    okscatter = okscatter && ( fabs(energy_hit3 - 1274.0) < 25.0 );
                }
            }

            // If Spatial Resolution:
            if (okscatter && doSmearHitPos)
            {
                CalculateSpatialSmearedHitPosition( PET_radius, hit1 );
                CalculateSpatialSmearedHitPosition( PET_radius, hit2 );
                CalculateSpatialSmearedHitPosition( PET_radius, hit3 );
            }

            // Do the 3-hit-TOF analysis
            if (okscatter)
            {
                DoThe3HitTOFStuff(iloop, hit1, hit2, hit3, lines, lines2nd, sourcePos
                    , scatterPosition1, scatterPosition2, scatterPosition3
                    , time_ps, gaussDist_time
                    , orig_hit1, orig_hit2, orig_hit3);
            }
        }           // if (distanceSrc2Hit2 > 0.0)  // if there is a 2nd hit on the LOR on the ring

        // new decay event time
        time_ps += delta_time_ps;

    }   // for (int iloop = 0; iloop < n_LORs; iloop++)  // loop over all events = LORs

    cout << "ngoodevents: " << ngoodevents << endl;
    cout << "ngoodAorBevents: " << ngoodAorBevents << endl;
    cout << "nogoodAnogoodB: " << nogoodAnogoodB << endl;

    gStyle->SetOptStat(1);
    TCanvas* m1 = new TCanvas("m1", " ", 700, 700);
    h1->SetMarkerStyle(20);
    h1->Draw("P");

    for (int iloop = 0; iloop < n_LORs; iloop++)
    {
        TLine* aLine = lines[iloop];
        if (aLine)
        {
            aLine->SetLineColor(kBlue);
            aLine->Draw("SAME");
        }

        TLine* aLine2 = lines2nd[iloop];
        if (aLine2)
        {
            aLine2->SetLineColor(kOrange);
            aLine2->Draw("SAME");
        }
    }
    m1->Print("hits.gif", "gif");

    TCanvas* m2 = new TCanvas("m2", " ", 1400, 700);
    gStyle->SetOptFit(1);   // moron!
        // Yes, ROOT, I want to see my fcking fit results in the histogram
        // Why the fck do you think I am doing the fit? Idiot!
    m2->Divide(2, 1);
    m2->cd(1);
    htsigma->Draw();
    htsigma->Fit("gaus", "Q");
    m2->cd(2);
    hdtsigma->Draw();
    hdtsigma->Fit("gaus", "Q");
    m2->Print("checking_generated_time_sigmas.gif", "gif");

    // -------------------------------------------------
    TCanvas* m3 = new TCanvas("m3", " ", 1400, 1400);
    gStyle->SetOptFit(1);   // moron!
    m3->Divide(2, 2);
    m3->cd(1);
    hErrorS_AB->Draw();
    hErrorS_AB->Fit("gaus", "Q");
    
    m3->cd(2);
    hErrorS_FA->Draw();
    hErrorS_FA->Fit("gaus", "Q");

    m3->cd(3);
    hErrorS_FB->Draw();
    hErrorS_FB->Fit("gaus", "Q");

    m3->cd(4);
    hErrorS_averaged->Draw();
    hErrorS_averaged->Fit("gaus", "Q");

    m3->Print("estimated_source_errors.gif", "gif");
    // -------------------------------------------------

    TCanvas* m4 = new TCanvas("m4", " ", 2100, 700);
    gStyle->SetOptFit(1);   // moron!
    m4->Divide(3, 1);
    m4->cd(1);
    hErrorS_AB_FA->Draw();
    hErrorS_AB_FA->Fit("gaus", "Q");

    m4->cd(2);
    hErrorS_AB_FB->Draw();
    hErrorS_AB_FB->Fit("gaus", "Q");

    m4->cd(3);
    hErrorS_FA_FB->Draw();
    hErrorS_FA_FB->Fit("gaus", "Q");

    m4->Print("estimated_source_errors_best2outof3.gif", "gif");

    // ugly ugly ugly
    AuxFinalizeH();

    return 0;
}

// ====================================================00

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

void CalculateHit2( const C3Vector& in_source, const C3Vector& in_direction, const double& in_lengthSrc2Hit2,
                    C3Vector& out_hit2 )
{
    out_hit2 = in_source + in_direction * (in_lengthSrc2Hit2/in_direction.GetLength());
}

// ====================================================00
// ====================================================00
// ====================================================00
// Here starts the 3 HITS TOF stuff....
//
void DoThe3HitTOFStuff(int iloop, const C3Vector& in_hit1, const C3Vector& in_hit2, const C3Vector& in_hit3
    , TLine** in_lines, TLine** in_lines2nd, const C3Vector& in_sourcePos
    , const C3Vector& in_scatterPosition1, const C3Vector& in_scatterPosition2, const C3Vector& in_scatterPosition3
    , const unsigned long long int& in_eventtime_ps, CVIPRandomGauss& in_gaussDist_time
    , const C3Vector& in_orig_hit1, const C3Vector& in_orig_hit2, const C3Vector& in_orig_hit3 )
{
    C3Vector hitA(in_hit1), hitB(in_hit2);
    // "A" is the hit closest to F (hit3)
    /* //
    {
        double distance_1_3 = (hit1 - hit3).GetLength();
        double distance_2_3 = (hit2 - hit3).GetLength();
        if (distance_2_3 < distance_1_3)
        // if (iloop%2 == 0)
        {
            hitA = hit2;
            hitB = hit1;
        }
    }
    */
    h1->Fill(hitA.GetX(), hitA.GetY());
    h1->Fill(hitB.GetX(), hitB.GetY());

    TLine* aLine = new TLine(hitA.GetX(), hitA.GetY(), hitB.GetX(), hitB.GetY());
    in_lines[iloop] = aLine;

    TLine* aLine2 = new TLine(in_sourcePos.GetX(), in_sourcePos.GetY(), in_hit3.GetX(), in_hit3.GetY());
    in_lines2nd[iloop] = aLine2;

    // Time calculation
    double dA = (hitA - in_sourcePos).GetLength();
    double dB = (hitB - in_sourcePos).GetLength();
    double d3 = (in_hit3 - in_sourcePos).GetLength();

    double dtA = dA/k_cspeed;
    double dtB = dB/k_cspeed;
    double dt3 = d3/k_cspeed;
    double timeA = in_eventtime_ps + dtA;
    double timeB = in_eventtime_ps + dtB;
    double time3 = in_eventtime_ps + dt3;

    double tmp_timeA( timeA );
    double tmp_timeB( timeB );
    double tmp_time3( time3 );

    ToPutInPlot toPutInPlot;

    toPutInPlot.time_A_clean = timeA;
    toPutInPlot.time_B_clean = timeB;
    toPutInPlot.time_F_clean = time3;

    toPutInPlot.dT_SA_clean = dtA;
    toPutInPlot.dT_SB_clean = dtB;
    toPutInPlot.dT_SF_clean = dt3;

    toPutInPlot.dT_FA_clean = time3 - timeA;
    toPutInPlot.dT_FB_clean = time3 - timeB;

    /*
    if (ngoodevents < k_plotGoodEvents)
    {
        cout << "CLEAN path -times: (1,2,3) " << dtA << " " << dtB << " " << dt3 << endl;
        cout << "CLEAN total times: " << timeA << " " << timeB << " " << time3 << endl;
    }
    */
    // TIME RESOLUTION
    // Smearing times !!!!!!!!!!
    //
    //  if (single_t_uncertainty_sigma_ps > 0)
    if (in_gaussDist_time.GetSigma() > 0)
    {
        double checkGauss1 = in_gaussDist_time.GetNewValue();
        timeA += checkGauss1;
        double checkGauss2 = in_gaussDist_time.GetNewValue();
        timeB += checkGauss2;
        double checkGauss3 = in_gaussDist_time.GetNewValue();
        time3 += checkGauss3;
    }

    toPutInPlot.time_A_smeared = timeA;
    toPutInPlot.time_B_smeared = timeB;
    toPutInPlot.time_F_smeared = time3;

    toPutInPlot.dT_SA_smeared = timeA - in_eventtime_ps;
    toPutInPlot.dT_SB_smeared = timeB - in_eventtime_ps;
    toPutInPlot.dT_SF_smeared = time3 - in_eventtime_ps;

    htsigma->Fill(timeA - tmp_timeA);
    htsigma->Fill(timeB - tmp_timeB);
    hdtsigma->Fill( (timeA-timeB) - (tmp_timeA-tmp_timeB));

    /*
    if (ngoodevents < k_plotGoodEvents)
    {
        cout << "SMEARED total times: " << timeA << " " << timeB << " " << time3 << endl;
    }
    */

    // The following source positions are based on Monte Carlo knowledge of travelling times
    C3Vector source_AF_travelling, source_BF_travelling;
    // The following source positions are estimations from using the time differences
    C3Vector source_AB_fromLOR, source_AF_iteration, source_BF_iteration;

    C3Vector lor_B2A = (hitA - hitB);
    //  cout << "S_AB <<< " << endl;
    {
        source_AB_fromLOR = CalculateSourcePosition( hitA, timeA, timeB, lor_B2A);

        // ERROR FROM S from A - B on LOR
        //
        double srcError = (source_AB_fromLOR - in_sourcePos).GetLength();

        // Which one is further away from hitA? (further along lor).
        // if        A --------- SOURCE ------------ S_AB    then error is negative
        // else if   A --------- S_AB -------------- SOURCE  then error is positive
        double dist1 = (source_AB_fromLOR - hitA).GetLength();
        double dist2 = (in_sourcePos - hitA).GetLength();
        int sign = (dist1 < dist2) ? 1 : -1;
        hErrorS_AB->Fill(srcError*sign);

        outputfile << " " << hitA << " " << timeA << " " << hitB << " " << timeB
                    << " " << in_hit3 << " " << time3
                    << " " << in_sourcePos << " " << source_AB_fromLOR << " " << srcError << endl;
    }
    //  cout << " >>> " << endl;

    // Moktar method (I) ***************************************************
    //
    /*
    C3Vector sourceS_F;
    bool isgood = MokhtarMethod_I(timeA, timeB, time3, hitA, hit3, lor, source_AB_fromLOR, toPutInPlot, sourceS_F);
    */

    // Moktar method (II) ***************************************************
    //
    C3Vector sourceS_FA;
    //  cout << " travel from A to B: " << endl;
    bool isgoodA = MokhtarTravelAlongLOR_II(lor_B2A, -1, hitA, in_hit3, time3 - timeA, sourceS_FA );
    C3Vector sourceS_FB;
    //  cout << " travel from B to A: " << endl;
    bool isgoodB = MokhtarTravelAlongLOR_II(lor_B2A,  1, hitB, in_hit3, time3 - timeB, sourceS_FB );

    if (!isgoodA && !isgoodB)
    {
        nogoodAnogoodB++;
    }

    toPutInPlot.dT_FA = time3 - timeA;
    toPutInPlot.dT_FB = time3 - timeB;
    toPutInPlot.length_AB = lor_B2A.GetLength();
    toPutInPlot.so_time_AB = lor_B2A.GetLength()/k_cspeed;
    toPutInPlot.Spos_AB = source_AB_fromLOR;
    toPutInPlot.Spos_FA = sourceS_FA;
    toPutInPlot.Spos_FB = sourceS_FB;

    if (isgoodA || isgoodB)
    {
        ngoodAorBevents++;
        if (isgoodA && isgoodB)
            ngoodevents++;

        // ERROR FROM S from FA
        //
        if (isgoodA) {
            double srcError = (sourceS_FA - in_sourcePos).GetLength();
            // The point closer to "hit A" is the most positive
            double dist1 = (sourceS_FA - hitA).GetLength();
            double dist2 = (in_sourcePos - hitA).GetLength();
            // if        A --------- SOURCE ------------ S_F     then error is negative
            // else if   A --------- S_F --------------- SOURCE  then error is positive

            int sign = (dist1 < dist2) ? 1 : -1;
            hErrorS_FA->Fill(srcError*sign);
        }
        // ERROR FROM S from FB
        //
        if (isgoodB) {
            double srcError = (sourceS_FB - in_sourcePos).GetLength();
            // The point closer to "hit A" is the most positive
            double dist1 = (sourceS_FB - hitA).GetLength();
            double dist2 = (in_sourcePos - hitA).GetLength();
            // if        A --------- SOURCE ------------ S_F     then error is negative
            // else if   A --------- S_F --------------- SOURCE  then error is positive

            int sign = (dist1 < dist2) ? 1 : -1;
            hErrorS_FB->Fill(srcError*sign);
        }

        // ERROR FROM AVERAGED S from LOR-AB AND F-A-B
        //
        {
            double divide_by_3(1.0/3.0);
            C3Vector sourceS_averaged;
            if (isgoodA && isgoodB)
                sourceS_averaged = (sourceS_FA + sourceS_FB + source_AB_fromLOR) * divide_by_3;
            else if (isgoodA)
                sourceS_averaged = (sourceS_FA + source_AB_fromLOR) * 0.5;
            else if (isgoodB)
                sourceS_averaged = (sourceS_FB + source_AB_fromLOR) * 0.5;
            else
            {
                cout << "HUH?" << isgoodA << " " << isgoodB << endl;
                exit(1);
            }
            double srcError = (sourceS_averaged - in_sourcePos).GetLength();
            double dist1 = (sourceS_averaged - hitA).GetLength();
            double dist2 = (in_sourcePos - hitA).GetLength();
            int sign = (dist1 < dist2) ? 1 : -1;
            hErrorS_averaged->Fill(srcError*sign);
        }
    }   // if (isgoodA || isgoodB)

    // if (isgoodA || isgoodB)
    {
        double D1 = (isgoodA) ? (sourceS_FA - source_AB_fromLOR).GetLength() : 100000;
        double D2 = (isgoodB) ? (sourceS_FB - source_AB_fromLOR).GetLength() : 100000;
        double D3 = (isgoodA && isgoodB) ? (sourceS_FB - sourceS_FA).GetLength() : 100000;
        {
            C3Vector sourceS_averaged;
            TH1D* histo(NULL);
            double srcError(0.0);
            if (D1 < D2 && D1 < D3)
            {
                sourceS_averaged = (sourceS_FA + source_AB_fromLOR) * 0.5;
                srcError = (sourceS_averaged - in_sourcePos).GetLength();
                histo = hErrorS_AB_FA;
            }
            else if (D2 < D1 && D2 < D3)
            {
                sourceS_averaged = (sourceS_FB + source_AB_fromLOR) * 0.5;
                srcError = (sourceS_averaged - in_sourcePos).GetLength();
                histo = hErrorS_AB_FB;
            }
            else if (D3 < D1 && D3 < D2)
            {
                sourceS_averaged = (sourceS_FB + sourceS_FA) * 0.5;
                srcError = (sourceS_averaged - in_sourcePos).GetLength();
                histo = hErrorS_FA_FB;
            }
            if (histo != NULL)
            {
                double dist1 = (sourceS_averaged - hitA).GetLength();
                double dist2 = (in_sourcePos - hitA).GetLength();
                int sign = (dist1 < dist2) ? 1 : -1;
                histo->Fill(srcError*sign);
            }
        }
    }   // no condition

    // if (iloop > 20)
    //  exit(1);

    // VISUALIZATION....
    //
    //  if (!isgoodA && !isgoodB && ngoodevents < k_plotGoodEvents)
    if (isgoodA && isgoodB && ngoodevents < k_plotGoodEvents)
    {
        gStyle->SetOptStat(0);
        TCanvas* m_temp = new TCanvas("m_temp", " ", 700, 700);

        h_petRing->Draw();

        petRing->SetLineWidth(2);
        petRing->SetLineColor(kBlack); // (kGreen+4);
        petRing->Draw("SAME");

        TLine* aLorLine = new TLine(hitA.GetX(), hitA.GetY(), hitB.GetX(), hitB.GetY());
        aLorLine->SetLineColor(kBlue);
        aLorLine->SetLineWidth(2);
        aLorLine->Draw("SAME");

        double radiusfactor(0.03);
		if ( PET_radius > 80.0 )
        	radiusfactor = 0.40;

        TEllipse* srcPoint = new TEllipse(in_sourcePos.GetX(), in_sourcePos.GetY(), radiusfactor*16, radiusfactor*16);
        srcPoint->SetLineColor(kGreen+2);
        srcPoint->SetFillColor(kGreen+2);
        srcPoint->Draw("SAME");

        TEllipse* Apoint = new TEllipse(hitA.GetX(), hitA.GetY(), radiusfactor*15, radiusfactor*15);
        Apoint->SetLineColor(kBlue+1);
        Apoint->SetFillColor(kBlue+1);
        Apoint->Draw("SAME");
        TText* text_A  = new TText(hitA.GetX(), hitA.GetY(), "A" );
        text_A->SetTextSize(0.03); text_A->SetTextColor(kBlue+1) ; text_A->Draw("SAME");

        TEllipse* Bpoint = new TEllipse(hitB.GetX(), hitB.GetY(), radiusfactor*15, radiusfactor*15);
        Bpoint->SetLineColor(kBlue);
        Bpoint->SetFillColor(kBlue);
        Bpoint->Draw("SAME");
        TText* text_B  = new TText(hitB.GetX(), hitB.GetY(), "B" );
        text_B->SetTextSize(0.03); text_B->SetTextColor(kBlue) ; text_B->Draw("SAME");

        TEllipse* Fpoint = new TEllipse(in_hit3.GetX(), in_hit3.GetY(), radiusfactor*15, radiusfactor*15);
        Fpoint->SetLineColor(kBlue-4);
        Fpoint->SetFillColor(kBlue-4);
        Fpoint->Draw("SAME");
        TText* text_F  = new TText(in_hit3.GetX(), in_hit3.GetY(), "F" );
        text_F->SetTextSize(0.03); text_F->SetTextColor(kBlue) ; text_F->Draw("SAME");

        TEllipse* Sabpoint = new TEllipse(source_AB_fromLOR.GetX(), source_AB_fromLOR.GetY(), radiusfactor*14, radiusfactor*14);
        Sabpoint->SetLineColor(kMagenta);
        Sabpoint->SetFillColor(kMagenta);
        TText* text_S_fromLOR  = new TText(source_AB_fromLOR.GetX(), source_AB_fromLOR.GetY()+0.4, "fromLOR" );
        text_S_fromLOR->SetTextSize(0.03); text_S_fromLOR->SetTextColor(kMagenta) ; text_S_fromLOR->Draw("SAME");
        Sabpoint->Draw("SAME");

        TEllipse* SFApoint = new TEllipse(sourceS_FA.GetX(), sourceS_FA.GetY(), radiusfactor*14, radiusfactor*14);
        SFApoint->SetLineColor(kRed);
        SFApoint->SetFillColor(kRed);
        TText* text_S_fromFA  = new TText(sourceS_FA.GetX(), sourceS_FA.GetY()-0.4, "fromFA" );
        text_S_fromFA->SetTextSize(0.03); text_S_fromFA->SetTextColor(kRed) ; text_S_fromFA->Draw("SAME");
        SFApoint->Draw("SAME");

        TEllipse* SFBpoint = new TEllipse(sourceS_FB.GetX(), sourceS_FB.GetY(), radiusfactor*14, radiusfactor*14);
        SFBpoint->SetLineColor(kPink+6);
        SFBpoint->SetFillColor(kPink+6);
        TText* text_S_fromFB  = new TText(sourceS_FB.GetX(), sourceS_FB.GetY()-1.2, "fromFB" );
        text_S_fromFB->SetTextSize(0.03); text_S_fromFB->SetTextColor(kPink+6) ; text_S_fromFB->Draw("SAME");
        SFBpoint->Draw("SAME");

		// All the scatter and smearing stuff
		//
        TEllipse* scatPoint1 = new TEllipse(in_scatterPosition1.GetX(), in_scatterPosition1.GetY()
            , radiusfactor*16, radiusfactor*16);
        scatPoint1->SetLineColor(kRed+2); scatPoint1->SetFillColor(kRed+2); scatPoint1->Draw("SAME");

        TEllipse* scatPoint2 = new TEllipse(in_scatterPosition2.GetX(), in_scatterPosition2.GetY()
            , radiusfactor*16, radiusfactor*16);
        scatPoint2->SetLineColor(kRed+2); scatPoint2->SetFillColor(kRed+2); scatPoint2->Draw("SAME");

        TEllipse* scatPoint3 = new TEllipse(in_scatterPosition3.GetX(), in_scatterPosition3.GetY()
            , radiusfactor*16, radiusfactor*16);
        scatPoint3->SetLineColor(kRed+2); scatPoint3->SetFillColor(kRed+2); scatPoint3->Draw("SAME");

        scatPoint1->Draw("SAME");
        scatPoint2->Draw("SAME");
        scatPoint3->Draw("SAME");

        TEllipse* Apoint_orig = new TEllipse(in_orig_hit1.GetX(), in_orig_hit1.GetY(), radiusfactor*15, radiusfactor*15);
        Apoint_orig->SetLineColor(kYellow); Apoint_orig->SetFillColor(kYellow+2);
        Apoint_orig->Draw("SAME");
        TText* text_A_orig  = new TText(in_orig_hit1.GetX(), in_orig_hit1.GetY(), "A" );
        text_A_orig->SetTextSize(0.03); text_A_orig->SetTextColor(kYellow+4) ; text_A_orig->Draw("SAME");

        TEllipse* Bpoint_orig = new TEllipse(in_orig_hit2.GetX(), in_orig_hit2.GetY(), radiusfactor*15, radiusfactor*15);
        Bpoint_orig->SetLineColor(kYellow); Bpoint_orig->SetFillColor(kYellow+2);
        Bpoint_orig->Draw("SAME");
        TText* text_B_orig  = new TText(in_orig_hit2.GetX(), in_orig_hit2.GetY(), "B" );
        text_B_orig->SetTextSize(0.03); text_B_orig->SetTextColor(kYellow+4) ; text_B_orig->Draw("SAME");

        TEllipse* Fpoint_orig = new TEllipse(in_orig_hit3.GetX(), in_orig_hit3.GetY(), radiusfactor*15, radiusfactor*15);
        Fpoint_orig->SetLineColor(kYellow); Fpoint_orig->SetFillColor(kYellow+2);
        Fpoint_orig->Draw("SAME");
        TText* text_F_orig  = new TText(in_orig_hit3.GetX(), in_orig_hit3.GetY(), "F" );
        text_F_orig->SetTextSize(0.03); text_F_orig->SetTextColor(kYellow+4) ; text_F_orig->Draw("SAME");

        TLine* aLorLineA = new TLine(hitA.GetX(), hitA.GetY(), in_sourcePos.GetX(), in_sourcePos.GetY());
        aLorLineA->SetLineColor(kBlue+2); aLorLineA->SetLineWidth(3);
        aLorLineA->Draw("SAME");
        
        TLine* aSrcToScat1 = new TLine(in_scatterPosition1.GetX(), in_scatterPosition1.GetY(), in_sourcePos.GetX(), in_sourcePos.GetY());
        aSrcToScat1->SetLineColor(kGreen+2); aSrcToScat1->SetLineWidth(3);
        aSrcToScat1->Draw("SAME");
        TLine* aScatToA = new TLine(in_scatterPosition1.GetX(), in_scatterPosition1.GetY(), hitA.GetX(), hitA.GetY());
        aScatToA->SetLineColor(kGreen+2); aScatToA->SetLineWidth(3);
        aScatToA->Draw("SAME");
        
        TLine* aSrcToScat2 = new TLine(in_scatterPosition2.GetX(), in_scatterPosition2.GetY(), in_sourcePos.GetX(), in_sourcePos.GetY());
        aSrcToScat2->SetLineColor(kGreen+2); aSrcToScat2->SetLineWidth(3);
        aSrcToScat2->Draw("SAME");
        TLine* aScatToB = new TLine(in_scatterPosition2.GetX(), in_scatterPosition2.GetY(), hitB.GetX(), hitB.GetY());
        aScatToB->SetLineColor(kGreen+2); aScatToB->SetLineWidth(3);
        aScatToB->Draw("SAME");

        TLine* aSrcToScat3 = new TLine(in_scatterPosition3.GetX(), in_scatterPosition3.GetY(), in_sourcePos.GetX(), in_sourcePos.GetY());
        aSrcToScat3->SetLineColor(kGreen+2); aSrcToScat3->SetLineWidth(3);
        aSrcToScat3->Draw("SAME");
        TLine* aScatToF = new TLine(in_scatterPosition3.GetX(), in_scatterPosition3.GetY(), in_hit3.GetX(), in_hit3.GetY());
        aScatToF->SetLineColor(kGreen+2); aScatToF->SetLineWidth(3);
        aScatToF->Draw("SAME");

        TString plot_filename = "drawing_event_NOtxt_";
        plot_filename += iloop;
        plot_filename += ".gif";
        m_temp->Print(plot_filename, "gif");

		// ADDITIONAL EXPLAINING TEXT
		//
        double top = PET_radius;
        double step = 0.05 * PET_radius;
        double xleft = -0.85 * PET_radius;

        TString stringDs = TString::Format("distances to real source (AS, BS, FS): %8.2f %8.2f %8.2f"
            , (hitA-in_sourcePos).GetLength(), (hitB-in_sourcePos).GetLength(), (in_hit3-in_sourcePos).GetLength());
        TText* infoDs = new TText(xleft, top, stringDs ); infoDs->SetTextSize(0.02); infoDs->Draw("SAME");

        TString stringT0 = TString::Format("clean times (A, B, F): %8.2f %8.2f %8.2f"
            , toPutInPlot.time_A_clean, toPutInPlot.time_B_clean, toPutInPlot.time_F_clean);
        TText* infoT0 = new TText(xleft,top-step, stringT0 ); infoT0->SetTextSize(0.02); infoT0->Draw("SAME");
        TString stringTS = TString::Format("measured times (A, B, F): %8.2f %8.2f %8.2f"
            , toPutInPlot.time_A_smeared, toPutInPlot.time_B_smeared, toPutInPlot.time_F_smeared);
        TText* infoTS = new TText(xleft,top-2*step, stringTS ); infoTS->SetTextSize(0.02); infoTS->Draw("SAME");

        TString stringM2 = TString::Format("dT SA clean: %4.2f", toPutInPlot.dT_SA_clean);
        TText* infoM2 = new TText(xleft,top-3*step, stringM2 ); infoM2->SetTextSize(0.02); infoM2->Draw("SAME");
        TString stringM1 = TString::Format("dT SB clean: %4.2f", toPutInPlot.dT_SB_clean);
        TText* infoM1 = new TText(xleft,top-4*step, stringM1 ); infoM1->SetTextSize(0.02); infoM1->Draw("SAME");
        TString string1 = TString::Format("dT SF clean: %4.2f", toPutInPlot.dT_SF_clean);
        TText* info1 = new TText(xleft, top-5*step, string1 ); info1->SetTextSize(0.02); info1->Draw("SAME");

        TString stringSAs = TString::Format("dT SA smeared: %4.2f", toPutInPlot.dT_SA_smeared);
        TText* infoSAs = new TText(xleft, top-6*step, stringSAs ); infoSAs->SetTextSize(0.02); 
            infoSAs->Draw("SAME");
        TString stringSBs = TString::Format("dT SB smeared: %4.2f", toPutInPlot.dT_SB_smeared);
        TText* infoSBs = new TText(xleft, top-7*step, stringSBs ); infoSBs->SetTextSize(0.02); 
            infoSBs->Draw("SAME");
        TString stringSFs = TString::Format("dT SF smeared: %4.2f", toPutInPlot.dT_SF_smeared);
        TText* infoSFs = new TText(xleft, top-8*step, stringSFs ); infoSFs->SetTextSize(0.02); 
            infoSFs->Draw("SAME");

        TString string2 = TString::Format("dT FA clean: %4.2f", toPutInPlot.dT_FA_clean);
        TText* info2 = new TText(xleft, top-9*step, string2 ); info2->SetTextSize(0.02); info2->Draw("SAME");
        TString string3 = TString::Format("dT FB clean: %4.2f", toPutInPlot.dT_FB_clean);
        TText* info3 = new TText(xleft, top-10*step, string3 ); info3->SetTextSize(0.02); info3->Draw("SAME");

        TString string5 = TString::Format("dT FA measured: %4.2f", toPutInPlot.dT_FA);
        TText* info5 = new TText(xleft, top-11*step, string5 ); info5->SetTextSize(0.02); info5->Draw("SAME");
        TString string6 = TString::Format("dT FB measured: %4.2f", toPutInPlot.dT_FB);
        TText* info6 = new TText(xleft, top-12*step, string6 ); info6->SetTextSize(0.02); info6->Draw("SAME");

        TString string7 = TString::Format("len AB: %4.2f", toPutInPlot.length_AB);
        TText* info7 = new TText(xleft, top-14*step, string7 ); info7->SetTextSize(0.02); info7->Draw("SAME");
        TString string8 = TString::Format("so time AB: %4.2f", toPutInPlot.so_time_AB);
        TText* info8 = new TText(xleft, top-15*step, string8 ); info8->SetTextSize(0.02); info8->Draw("SAME");

        double len_SA = (hitA - sourceS_FA).GetLength();
        double time_SA = len_SA/k_cspeed;
        TString stringSA = TString::Format("iteration SA len: %4.2f time: %4.2f", len_SA, time_SA);
        TText* infoSA = new TText(xleft, top-16*step, stringSA ); infoSA->SetTextSize(0.02); 
            infoSA->Draw("SAME");

        double len_SFA = (in_hit3 - sourceS_FA).GetLength();
        double time_SFA = len_SFA/k_cspeed;
        TString stringSFA = TString::Format("iteration(A) SF len: %4.2f time: %4.2f", len_SFA, time_SFA);
        TText* infoSFA = new TText(xleft, top-17*step, stringSFA ); infoSFA->SetTextSize(0.02); 
            infoSFA->Draw("SAME");

        double iter_dT_F_A = time_SFA - time_SA;
        double errorA = fabs(iter_dT_F_A - toPutInPlot.dT_FA);
        TString string_iter_dtFA = TString::Format("iteration dT FA: %4.2f error: %4.2f", iter_dT_F_A, errorA);
        TText* info_iter_dtFA = new TText(xleft, top-18*step, string_iter_dtFA ); 
            info_iter_dtFA->SetTextSize(0.02);
            info_iter_dtFA->Draw("SAME");

        double len_SB = (hitB - sourceS_FB).GetLength();
        double time_SB = len_SB/k_cspeed;
        TString stringSB = TString::Format("iteration SB len: %4.2f time: %4.2f", len_SB, time_SB);
        TText* infoSB = new TText(xleft,  top-19*step, stringSB ); infoSB->SetTextSize(0.02); 
            infoSB->Draw("SAME");

        double len_SFB = (in_hit3 - sourceS_FB).GetLength();
        double time_SFB = len_SFB/k_cspeed;
        TString stringSFB = TString::Format("iteration(B) SF len: %4.2f time: %4.2f", len_SFB, time_SFB);
        TText* infoSFB = new TText(xleft, top-20*step, stringSFB ); infoSFB->SetTextSize(0.02); 
            infoSFB->Draw("SAME");

        double iter_dT_F_B = time_SFB - time_SB;
        double errorB = fabs(iter_dT_F_B - toPutInPlot.dT_FB);
        TString string_iter_dtFB = TString::Format("iteration dT FB: %4.2f error: %4.2f", iter_dT_F_B, errorB);
        TText* info_iter_dtFB = new TText(xleft, top-21*step, string_iter_dtFB ); 
            info_iter_dtFB->SetTextSize(0.02);
            info_iter_dtFB->Draw("SAME");

        TString string13 = TString::Format("S from AB: %4.2f %4.2f", toPutInPlot.Spos_AB.GetX(), toPutInPlot.Spos_AB.GetY());
        TText* info13 = new TText(xleft, top-22*step, string13 ); info13->SetTextSize(0.02); 
            info13->Draw("SAME");
        TString string14 = TString::Format("S from FA: %4.2f %4.2f", toPutInPlot.Spos_FA.GetX(), toPutInPlot.Spos_FA.GetY());
        TText* info14 = new TText(xleft, top-23*step, string14 ); info14->SetTextSize(0.02); 
            info14->Draw("SAME");
        TString string15 = TString::Format("S from FB: %4.2f %4.2f", toPutInPlot.Spos_FB.GetX(), toPutInPlot.Spos_FB.GetY());
        TText* info15 = new TText(xleft, top-24*step, string15 ); info15->SetTextSize(0.02); 
            info15->Draw("SAME");
        TString string16 = TString::Format("S real: %4.2f %4.2f", in_sourcePos.GetX(), in_sourcePos.GetY());
        TText* info16 = new TText(xleft, top-25*step, string16 ); info16->SetTextSize(0.02); 
            info16->Draw("SAME");

        // MORE CRAP
        top = -20.0;
        xleft = 0.0;

        TString posAStr = TString::Format("posA: %4.2f %4.2f", hitA.GetX(), hitA.GetY());
        TText* posAtxt = new TText(xleft, top, posAStr ); posAtxt->SetTextSize(0.02); posAtxt->Draw("SAME");
        TString posBStr = TString::Format("posB: %4.2f %4.2f", hitB.GetX(), hitB.GetY());
        TText* posBtxt = new TText(xleft, top-step, posBStr ); posBtxt->SetTextSize(0.02); posBtxt->Draw("SAME");
        TString posFStr = TString::Format("posF: %4.2f %4.2f", in_hit3.GetX(), in_hit3.GetY());
        TText* posFtxt = new TText(xleft, top-2*step, posFStr ); posFtxt->SetTextSize(0.02);
            posFtxt->Draw("SAME");

        // ---
        TString posSLORStr = TString::Format("posSLOR: %4.2f %4.2f", source_AB_fromLOR.GetX(), source_AB_fromLOR.GetY());
        TText* posSLORtxt = new TText(xleft, top-3*step, posSLORStr ); posSLORtxt->SetTextSize(0.02);
            posSLORtxt->Draw("SAME");

        double fckingPieceOfShit_distance_S_LOR_TO_A = (source_AB_fromLOR - hitA).GetLength();
        double so_tmp_SLOR_to_A = fckingPieceOfShit_distance_S_LOR_TO_A/k_cspeed;
        TString fckA = TString::Format("Distance SLOR TO A: %4.2f time: %4.2f", fckingPieceOfShit_distance_S_LOR_TO_A, so_tmp_SLOR_to_A);
        TText* fckAtxt = new TText(xleft, top-4*step, fckA ); fckAtxt->SetTextSize(0.02);
            fckAtxt->Draw("SAME");

        double fckingPieceOfShit_distance_S_LOR_TO_B = (source_AB_fromLOR - hitB).GetLength();
        double so_tmp_SLOR_to_B = fckingPieceOfShit_distance_S_LOR_TO_B/k_cspeed;
        TString fckB = TString::Format("Distance SLOR TO B: %4.2f time: %4.2f", fckingPieceOfShit_distance_S_LOR_TO_B, so_tmp_SLOR_to_B);
        TText* fckBtxt = new TText(xleft, top-5*step, fckB ); fckBtxt->SetTextSize(0.02);
            fckBtxt->Draw("SAME");
            
        // ---
        TString posS_FAStr = TString::Format("posS_from_dt_FA: %4.2f %4.2f", sourceS_FA.GetX(), sourceS_FA.GetY());
        TText* posS_FAtxt = new TText(xleft, top-6*step, posS_FAStr ); posS_FAtxt->SetTextSize(0.02); 
            posS_FAtxt->Draw("SAME");

        double fckingPieceOfShit_distance_S_FA_TO_A = (sourceS_FA - hitA).GetLength();
        double so_tmp_SFA_to_A = fckingPieceOfShit_distance_S_FA_TO_A/k_cspeed;
        TString fck_FA_A = TString::Format("Distance S_FA TO A: %4.2f time: %4.2f", fckingPieceOfShit_distance_S_FA_TO_A, so_tmp_SFA_to_A);
        TText* fck_FA_Atxt = new TText(xleft, top-7*step, fck_FA_A ); fck_FA_Atxt->SetTextSize(0.02); 
            fck_FA_Atxt->Draw("SAME");

        plot_filename = "drawing_event_txt_";
        plot_filename += iloop;
        plot_filename += ".gif";
        m_temp->Print(plot_filename, "gif");

        delete m_temp;
        m_temp = 0;

        // ======================= =================== ========
    }
}


// ------------------------------------------------------------------------------
// -- new stuff 2019-09-19

double GetSingleTimeResolutionGauss()
{
    // Time resolution per single time, instead of time-difference
    // (t_a - t_b) +- Gauss_FWHM( time12_uncertainty_ps)
    // d_t^2 + d_t^2 = sigma^2
    // 2 d_t^2 = sigma^2
    // d_t^2 = 0.5 * sigma^2
    // d_t = sqrt(0.5) * sigma

    double d_t = 0.5 * sqrt(2) * time12_uncertainty_ps_sigma;
    return d_t;
}

// -------------------------------------------------------------

// Hit2 might not be lying on the LOR
// S_AB = M - c * dT_AB/2 * unit_vector_along_LOR
// M = (A + B)/2 = B + LOR/2 = A - LOR/2

C3Vector CalculateSourcePosition( const C3Vector& in_hit1, const double& in_time_1,
                                  const double& in_time_2, const C3Vector& in_lor_from2_to_1 )
{
    bool debug(false);

    double delta_time = (in_time_1 - in_time_2);
    double del_distance_from_time = k_cspeed * delta_time / 2.0;

    C3Vector lor_direction(in_lor_from2_to_1);
    lor_direction.Normalize(1.0);
    C3Vector D_along_lor = lor_direction * del_distance_from_time;
    if (debug) cout << "    delta_time: " << delta_time << " ; Distance D from time: " << del_distance_from_time
                    << " ; D_along_lor: " <<  D_along_lor << endl;

    C3Vector center = in_hit1 - in_lor_from2_to_1 * 0.5;
    C3Vector sourceGuess = center - D_along_lor;
    if (debug) cout << "    center: " << center << " sourceGuess: " << sourceGuess << endl;

    return sourceGuess;
}

// -----------------------------------------

double GetCosAlphaThreeTriangleSides( const double& in_side1_AF, const double& in_side2_AS, const double& in_side3_FS )
{
    // in_side3_FS^2 = in_side1_AF^2 + in_side2_AS^2 - 2 * in_side1_AF * in_side2_AS * cos(alpha)
    // 2 * in_side1_AF * in_side2_AS * cos(alpha) = in_side1_AF^2 + in_side2_AS^2 - in_side3_FS^2
    // cos(alpha) = cos(alpha_12) = (in_side1_AF^2 + in_side2_AS^2 - in_side3_FS^2) / (2 * in_side1_AF * in_side2_AS)

    double cos_alpha = (in_side1_AF*in_side1_AF + in_side2_AS*in_side2_AS - in_side3_FS*in_side3_FS);
    cos_alpha = cos_alpha/ (2 * in_side1_AF * in_side2_AS);
    return cos_alpha;
}

// --------------------------------_

double GetDeltaTimeDifference(const C3Vector& approx_source, const C3Vector& in_hit1, const double& in_time1
    , const C3Vector& in_hit3, const double& in_time3)
{
    C3Vector AS = in_hit1 - approx_source;
    C3Vector FS = in_hit3 - approx_source;
    double delta_time_measured = (in_time1 - in_time3);
    double delta_time_from_path = AS.GetLength()/k_cspeed - FS.GetLength()/k_cspeed;

cout << "AS: " << AS << " path_AS: " << AS.GetLength() << " ~time_AS: " << AS.GetLength()/k_cspeed << endl;
cout << "FS: " << FS << " path_FS: " << FS.GetLength() << " ~time_FS: " << FS.GetLength()/k_cspeed << endl;
cout << "delta_time_measured: " << delta_time_measured << " delta_time_from_path: " << delta_time_from_path << endl;

    double factor = (delta_time_measured - delta_time_from_path);
    return factor;
}

// ----------------------------------

double UpdateApproximateSource(C3Vector& approx_source, const C3Vector& in_hit1, const double& in_time1
    , const C3Vector& in_hit3, const double& in_time3, const C3Vector& in_lor )
{
    bool debug(true);
    if (debug) cout << "UPDATING: " << endl;

    double previous_error = GetDeltaTimeDifference(approx_source, in_hit1, in_time1, in_hit3, in_time3);

    if (debug) cout << "previous_error (time): " << previous_error << endl;
    double factor = previous_error * k_cspeed;
    // cout << "factor (space): " << factor << endl;
    C3Vector delta_S(in_lor);
    delta_S.Normalize(factor);

    approx_source = approx_source - delta_S;

    // cout << " factor: " << factor << " delta_S: " << delta_S
    //      << " new approx S: " << approx_source << endl;

    if (debug) cout << " NEW ITERATED SOURCE: " << approx_source << endl;

    return previous_error;
}

// ------------------------

void
GetMonteCarloKnowledgeSourcePositionFromTravellingTimes( C3Vector& out_sourceTravelling, const double& in_time_event_ps,
    const C3Vector& in_hit_LOR, const double& in_time_LORhit, const C3Vector& in_hit_3, const double& in_time_3 )
{
    double d1_travelling = in_time_LORhit - in_time_event_ps;
    double d3_travelling = in_time_3 - in_time_event_ps;
    double path_AS_travelling = d1_travelling * k_cspeed;
    double path_FS_travelling = d3_travelling * k_cspeed;
    double path_AF_real = (in_hit_LOR - in_hit_3).GetLength();

    double cos_alpha_FaS = GetCosAlphaThreeTriangleSides( path_AF_real, path_AS_travelling, path_FS_travelling );

    // Even though it's in 3 dimensions, for the moment, we only use 2 dimensions.
    //
    double alpha_FaS = std::acos(cos_alpha_FaS);
    // cout << "cos_alpha_FaS: " << cos_alpha_FaS << " alpha (degrees): " << alpha_FaS*180.0/kPI << endl;
    C3Matrix rotationMatrix( Axis_Z, -1.0*alpha_FaS);
    C3Vector aVector = in_hit_3 - in_hit_LOR;         // vector from A to F
    aVector = aVector * rotationMatrix;     // projected unto AS
    aVector.Normalize( path_AS_travelling );     // normalized to length of AS
    out_sourceTravelling = in_hit_LOR + aVector;       // Position of "source" of A-S-F triangle

    bool debug(true);
    if (debug)
    {
        cout << "Monte Carlo knowledge: <<<" << endl;
        cout << " path_AS(based on known dtime): " << path_AS_travelling
             << " ~time_AS: " << path_AS_travelling/k_cspeed << " check (known travel time): " << d1_travelling << endl;
        cout << " path_FS(based on known dtime): " << path_FS_travelling
             << " ~time_FS: " << path_FS_travelling/k_cspeed << " check (known travel time): " << d3_travelling << endl;
        double delta_time_from_path = path_AS_travelling/k_cspeed - path_FS_travelling/k_cspeed;
        cout << " Delta t from Monte Carlo known travelling path: " << delta_time_from_path << endl;
        cout << "S_AF position: " << out_sourceTravelling << endl;
        cout << ">>> MC knowledge" << endl;
    }
}

// --------------------------------------------------

void GetIteratedSourceApproximation( C3Vector& out_source_XF_iteration, C3Vector& in_other_hit_fromLOR
    , const C3Vector& in_source_XF_travelling, const C3Vector& in_realSourcePos
    , const C3Vector& in_hit_LOR, const double& in_time_LORhit, const C3Vector& in_hit_3, const double& in_time_3
    , TH1D* h_S_XF_srcDeviation, TH1D* h_XsF_timeDeviation, TH1D* hErrorS_AX )
{
    // Approximate
    //

    C3Vector lor = (in_hit_LOR - in_other_hit_fromLOR);

    //  source_XF_iteration = source_AB_fromLOR; // initial guess
    out_source_XF_iteration = in_other_hit_fromLOR; // initial guess
    cout << "START Source_AF: " << out_source_XF_iteration << endl;
    double err(0.0), srcDev(0.0);
    int i(0);

    const int NITER(1);
    for (i = 0; i < NITER; i++)
    {
        cout << i << endl;
        srcDev = (out_source_XF_iteration - in_source_XF_travelling).GetLength();
        h_S_XF_srcDeviation->Fill(i, srcDev);
        err = UpdateApproximateSource( out_source_XF_iteration, in_hit_LOR, in_time_LORhit, in_hit_3, in_time_3, lor );
        h_XsF_timeDeviation->Fill(i, err);
    }
    srcDev = (out_source_XF_iteration - in_source_XF_travelling).GetLength();
    h_S_XF_srcDeviation->Fill(i, srcDev);
    err = GetDeltaTimeDifference(out_source_XF_iteration, in_hit_LOR, in_time_LORhit, in_hit_3, in_time_3);
    h_XsF_timeDeviation->Fill(i, err);

    double srcError = (out_source_XF_iteration - in_realSourcePos).GetLength();
    // cout << std::setprecision (10)
    if (srcError < 1e-5) srcError = 0.0;   // SOMETIMES C++ is a bloody moron
    int sign = (out_source_XF_iteration.GetX() > in_realSourcePos.GetX()) ? -1 : 1;
    hErrorS_AX->Fill(srcError*sign);

}

// =======================

bool MokhtarTravelAlongLOR_II(const C3Vector& in_lor_B2A, int in_direction, const C3Vector& in_lorHit
    , const C3Vector& in_hit3, const double& in_delta_T_FX_measured
    , C3Vector& io_sourceS_F )
{
    C3Vector lor_B2A_unit(in_lor_B2A);
    lor_B2A_unit.Normalize();

    // Going from A to B or from B to A
    double segment = 0.1;
    double goodsegment = 0.0;
    double error(0);
    double preverror = 10000;
    C3Vector prev_S;
    //  bool stop(false);
    while (segment < in_lor_B2A.GetLength()) // && !stop)
    {
        io_sourceS_F = in_lorHit + ( lor_B2A_unit * in_direction * segment);

        double len_FS = (io_sourceS_F - in_hit3).GetLength();
        double delta_T_XS = segment / k_cspeed;
        double delta_T_FS = len_FS / k_cspeed;
        double guess_delta_T_FX = delta_T_FS - delta_T_XS;

        double error = fabs(guess_delta_T_FX - in_delta_T_FX_measured);
        /*
        if (error < 2)
        {
            cout << " segment: " << segment << " FS: " << len_FS << " guess T_Fx: " << guess_delta_T_FX
                 << " S: " <<  io_sourceS_F << " error: " << error << endl;
        }
        */
        if (error < preverror)
        {
            prev_S = io_sourceS_F;
            preverror = error;
            goodsegment = segment;
        }
        // else
        //    stop = true;
        segment += 0.1;
    }

    // bool isgood = stop && segment > 2.0;
        // if stop == true, it means we found a minimal error
    bool isgood = goodsegment > 1.0 && goodsegment < (in_lor_B2A.GetLength() - 1.0);
    //  if (isgood)
    {
        if (isgood)
        {
            io_sourceS_F = prev_S;
            //  cout << "final source: " << io_sourceS_F << endl;
        }
        else
        {
            io_sourceS_F = in_lorHit;
        }
        //  cout << "Final error: " << preverror << endl;
    }

    return isgood;
    //  return true;
}

// =========================

bool MokhtarMethod_I(const double& in_timeA, const double& in_timeB, const double& in_time3
    , const C3Vector& in_hitA, const C3Vector& in_hit3, const C3Vector& in_lor, const C3Vector& in_source_AB_fromLOR
    , ToPutInPlot& io_toPutInPlot, C3Vector& io_sourceS_F)
{
    bool isgood(false);
    double solution1, solution2;
    double distance_SA, distance_SA_2;
    {
        double delta_t_FA = in_time3 - in_timeA;
        double delta_t_FB = in_time3 - in_timeB;

        double lorLength = in_lor.GetLength();
        double tSA_plus_t_SB = lorLength/k_cspeed;
        double t_SF = (delta_t_FA + delta_t_FB + tSA_plus_t_SB)/2.0;
        double distance_SF = k_cspeed * t_SF;

        double distance_AF = (in_hitA - in_hit3).GetLength();
        double cosAngle_SF = CalculateCosAngleBD( in_lor, (in_hit3 - in_hitA) );
            // the angle opposite side SF = the angle between AS and AF = angle between lor and AF
        int numsolutions = CalculateLengthBOfTriangle(distance_SF, distance_AF, cosAngle_SF, solution1, solution2 );
            /*
            In ANY triangle with sides with lengths SF, AS, AF,
            and an angle alpha_opposite_SF between AS(unknown, but along known lor) and AF,
            the following is true:
                AS^2 - 2 * AF * cos(alpha_opposite_SF)* AS + (AF^2 - SD^2) = 0
            Quadratic equation:
                unknown AS -> x
                    a = 1
                    b = - 2 * C * cos(alpha_opposite_A)
                    c = (C^2 - A^2)
            */

        distance_SA = 0.0;
        distance_SA_2 = -1.0;
        double diff1(0.0), diff2(0.0);
        C3Vector tmpV1, tmpV2;
        if (numsolutions == 1)
            distance_SA = solution1;
        else if (numsolutions == 2)
        {
            C3Vector lor1(in_lor);
            lor1.Normalize(solution1);
            tmpV1 = in_hitA + lor1;

            C3Vector lor2(in_lor);
            lor2.Normalize(solution2);
            tmpV2 = in_hitA + lor2;

            diff1 = (tmpV1 - in_source_AB_fromLOR).GetLength();
            diff2 = (tmpV2 - in_source_AB_fromLOR).GetLength();
            if (diff1 <= diff2)
            {
                distance_SA = solution1;
                distance_SA_2 = solution2;
                io_sourceS_F = tmpV1;
            }
            else
            {
                distance_SA = solution2;
                distance_SA_2 = solution1;
                io_sourceS_F = tmpV2;
            }
        }
        /*
        if (ngoodevents < k_plotGoodEvents)
        {
            cout << iloop << " # solutions: " << numsolutions << " (1): " << solution1 << " (2): " << solution2
                << " dSA: " << distance_SA << " (2): " << distance_SA_2 << endl;
            cout << "because, diff1: " << diff1 << " diff2: " << diff2
                << " coming from solvec1: " << tmpV1 << " or 2: " << tmpV2 << " and S_ab: " << source_AB_fromLOR << endl;
        }
        */
        if (numsolutions > 0)
        {
            isgood = true;

            io_toPutInPlot.dT_FA = delta_t_FA;
            io_toPutInPlot.dT_FB = delta_t_FB;
            io_toPutInPlot.length_AB = lorLength;
            io_toPutInPlot.so_time_AB = tSA_plus_t_SB;
            io_toPutInPlot.calculated_time_SF = t_SF;
            io_toPutInPlot.so_length_SF = distance_SF;
            // io_toPutInPlot.calculated_length_AS_1 = distance_SA;
            // io_toPutInPlot.calculated_length_AS_2 = distance_SA_2;

            io_toPutInPlot.Spos_AB = in_source_AB_fromLOR;
            io_toPutInPlot.Spos_FA = io_sourceS_F;
        }
    }

    return isgood;
}


// -----------------------------------------
// https://www.researchgate.net/publication/272440212_Teaching_Mathematics_with_Mathematical_Software/figures?lo=1
void TestFormulas()
{
    double R = 28;
    double S = 12;
    double cos_alpha = -0.5;

    double sol1(-1), sol2(-1);

    // CalculateDistanceSrcToHit2(R, S, cos_alpha, sol1, sol2);
    CalculateLengthBOfTriangle(R, S, cos_alpha, sol1, sol2);
    cout << "solutions: " << sol1 << " " << sol2 << endl;
    cout << "R = a: " << R << " S = b: " << S << " sol1 = c: " << sol1 << " cos_a: " << cos_alpha << endl;

    // ===========================
    double cos_alpha_00 = GetCosAlphaThreeTriangleSides( S, sol1, R );
    cout << "cos_alpha_00: " << cos_alpha_00 << endl;

    // ===========================
    {
        cout << endl;
        double side_1 = 3;
        double side_2 = 4;
        double side_3 = 5;
        double cos_a_12 = GetCosAlphaThreeTriangleSides( side_1, side_2, side_3 );
        cout << "cos_a_12: " << cos_a_12 << " degrees: " << acos(cos_a_12) * 180.0/kPI << endl;

        double cos_a_13 = GetCosAlphaThreeTriangleSides( side_1, side_3, side_2 );
        cout << "cos_a_13: " << cos_a_13 << " degrees: " << acos(cos_a_13) * 180.0/kPI << endl;

        double cos_a_23 = GetCosAlphaThreeTriangleSides( side_2, side_3, side_1 );
        cout << "cos_a_23: " << cos_a_23 << " degrees: " << acos(cos_a_23) * 180.0/kPI << endl;
    }
    {
        cout << endl;
        double side_1 = 20;
        double side_2 = 12;
        double side_3 = 28;
        double cos_a_12 = GetCosAlphaThreeTriangleSides( side_1, side_2, side_3 );
        cout << "cos_a_12: " << cos_a_12 << " degrees: " << acos(cos_a_12) * 180.0/kPI << endl;

        double cos_a_13 = GetCosAlphaThreeTriangleSides( side_1, side_3, side_2 );
        cout << "cos_a_13: " << cos_a_13 << " degrees: " << acos(cos_a_13) * 180.0/kPI << endl;

        double cos_a_23 = GetCosAlphaThreeTriangleSides( side_2, side_3, side_1 );
        cout << "cos_a_23: " << cos_a_23 << " degrees: " << acos(cos_a_23) * 180.0/kPI << endl;
    }
    {
        cout << endl;
        double side_1 = 180;
        double side_2 = 238;
        double side_3 = 340;
        double cos_a_12 = GetCosAlphaThreeTriangleSides( side_1, side_2, side_3 );
        cout << "cos_a_12: " << cos_a_12 << " degrees: " << acos(cos_a_12) * 180.0/kPI << endl;

        double cos_a_13 = GetCosAlphaThreeTriangleSides( side_1, side_3, side_2 );
        cout << "cos_a_13: " << cos_a_13 << " degrees: " << acos(cos_a_13) * 180.0/kPI << endl;

        double cos_a_23 = GetCosAlphaThreeTriangleSides( side_2, side_3, side_1 );
        cout << "cos_a_23: " << cos_a_23 << " degrees: " << acos(cos_a_23) * 180.0/kPI << endl;
    }
    cout << endl;
}























