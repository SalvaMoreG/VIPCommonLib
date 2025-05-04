#include "CVIPToyPET_aux.h"

#include <iostream>
#include <cmath>

#include "CVIPUtils.h"
#include "CVIPRandom.h"
#include "CVIP3Vector.h"
#include "CVIP3Matrix.h"
#include "CVIPPhysicsUtils.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

const double P_scatter(1.0);    // between 0 and 1.0

TH2D* h_EvsAngle(0);

bool CalculateSpatialSmearedHitPosition( const double& in_radiusPET, C3Vector& io_hitPosition )
{
    const double voxelsize = 2.0; // mm
    int number_of_sections = int(2.0 * kPI * in_radiusPET/voxelsize);
    // cout << "#sections in PET ring: " << number_of_sections << endl;
    double delta_angle = 2.0 * kPI/number_of_sections;
    //  cout << "delta_angle: " << delta_angle << endl;

    double angle(0.0);
    bool ok = VIPUtils::GetPositionAngle2D( io_hitPosition.GetX(), io_hitPosition.GetY(), angle);
    if (!ok) return false;

    // cout << "angle hit location (BEFORE): " << angle << " " << angle * 180/kPI
    //  << " position: " << io_hitPosition << endl;
    int section_number = int(angle/delta_angle);
    //  cout << "So, the current section number: " << section_number << endl;
    angle = ((float)section_number + 0.5) * delta_angle;
    //  cout << "So, new angle: " << angle << endl;

    io_hitPosition.SetX( in_radiusPET * std::cos(angle));
    io_hitPosition.SetY( in_radiusPET * std::sin(angle));
    //  cout << " NEW position: " << io_hitPosition << endl;

	return true;
}

// ===========================================================================

bool CalculateScatteredTrack( const double& in_radiusPET, const double& in_radiusSphere
           , const C3Vector& in_sourcePosition, C3Vector& io_hitPosition, C3Vector& out_scatterPosition
           , double& io_newEnergy )
{
    out_scatterPosition = in_sourcePosition;

    if (h_EvsAngle == 0)
    {
        h_EvsAngle = new TH2D("h_EvsAngle", "E vs angle", 50, 0.0, 50.0, 1275, 0.0, 1275.0);
//        h_EvsAngle = new TH2D("h_EvsAngle", "E vs angle", 50, 0.0, 50.0, 100, 450.0, 550.0);
    }

    // 1.0 Calculate the probability that the track will scatter
    CVIPRandomUniform uniformDist(0.0, 1.0);
    double h_or_m = uniformDist.GetNewValue();
    if (h_or_m > P_scatter)
        return false;

    // 2.0 Calculate where in the sphere the gamma scatters
    double distance = uniformDist.GetNewValue();
    distance = distance * in_radiusSphere;
distance = 2.0;
    C3Vector old_direction = io_hitPosition - in_sourcePosition;
    old_direction.Normalize();
    out_scatterPosition = in_sourcePosition + old_direction * distance;

    // 3.0 Using Klein-Nishina calculate the scattering angle
    bool accepted = false;
    double scatterTheta(0.0);
//    /*
    double maxAngle = 15.0 * kPI/180.0;
    while (!accepted)
    {
        scatterTheta = uniformDist.GetNewValue() * maxAngle;
        double Pscatter = CPhysicsUtils::CalculateKleinNishinaProb( mass_electron_keV, cos(scatterTheta) );

        h_or_m = uniformDist.GetNewValue();
        if (h_or_m < Pscatter)
            accepted = true;
    }
//    */
//  scatterTheta = kPI * 10.0/180.0;
    //  cout << "scatterTheta: " << scatterTheta << " degrees: " << 180.0 * scatterTheta/kPI << endl;

    // 4. Calculate the direction of the vector from scatter position to new hit B'
    // i.e. rotate the vector with direction B - S', with scatterTheta, to direction B' - S'
    C3Matrix rotMatrix( Axis_Z, scatterTheta);
    C3Vector new_direction = old_direction * rotMatrix;
    //  cout << "old_direction: " << old_direction << " new: " << new_direction << endl;

    // 5. Calculate the angle between new direction and the line from S' to the PET-center
    C3Vector vec_Scat_to_O = C3Vector(0, 0, 0) - out_scatterPosition;
    double cos_alpha = new_direction.GetScalarProductCosAngle( vec_Scat_to_O );
    //  cout << "cos_alpha: " << cos_alpha << " alpha: " << acos(cos_alpha)
    //       << " degrees: " << 180.0 * acos(cos_alpha)/kPI << endl;

    // 6. Calculate the length of the new direction, from S' until unknown B'
    // Using this function:
    // triangle with sides A = R_sphere = (B' - O); B = (B' - S') and C (S' - O)
    // unknown: length B' - S'
    double scatterpos_to_origin_len = out_scatterPosition.GetLength();
        // the origin is at (0, 0, 0), so no need to subtract
    double new_length1, new_length2;
    //  cout << "INPUT CalcBLen, R: " << in_radiusPET << " (S'-O): " << scatterpos_to_origin_len
    //       << " cos_alpha: " << cos_alpha << endl;
    int number_of_solutions = CalculateLengthBOfTriangle(in_radiusPET, scatterpos_to_origin_len
                                    , cos_alpha, new_length1, new_length2);

    if (number_of_solutions == 0)
        return false;
    else if (number_of_solutions == 1)
        new_direction = new_direction * new_length1;
    else if (new_length1 > 0 && new_length2 < 0)
        new_direction = new_direction * new_length1;
    else if (new_length1 < 0 && new_length2 > 0)
        new_direction = new_direction * new_length2;
    else if (new_length1 < 0 && new_length2 < 0)
    {
        cout << number_of_solutions << " " << new_length1 << " " << new_length2 << endl;
        new_direction = new_direction * new_length2;
    }
    else if (new_length1 < 0 && new_length2 < 0)
        return false;

    // 7 Calculate new hit B' position
    io_hitPosition = out_scatterPosition + new_direction;

    // 8. Also calculate new energy of scattered gamma

    double totEnergy( io_newEnergy );
    bool ok = VIPUtils::GetComptonEnergyFromAngleRad(scatterTheta, totEnergy, io_newEnergy );
        // returns E1 (energy lost in the scattering)

    // Get E2
    io_newEnergy = totEnergy - io_newEnergy;
    //  cout << "final out_energy: " << out_newEnergy << endl;

    //  cout << "filling h_EvsAngle, angle: " << scatterTheta*180.0/kPI << " E': " << out_newEnergy << endl;
    h_EvsAngle->Fill(scatterTheta*180.0/kPI, io_newEnergy);

    return true;
}

void AuxFinalizeH()
{
    if (h_EvsAngle)
    {
        TCanvas* tmpaux = new TCanvas("tmpaux", " ", 700, 700);
        h_EvsAngle->SetMarkerStyle(20);
        h_EvsAngle->Draw("P");
        tmpaux->Print("EvsAngle.gif", "gif");
    }
}

// ****************************************************************
// ****************************************************************
// Explanation of how to get two hits on a PET ring, given a source which is NOT at the centre (0, 0, 0)
// See:
// https://www.researchgate.net/publication/272440212_Teaching_Mathematics_with_Mathematical_Software/figures?lo=1
// doc/Triangle_Solution_Quadratic.png
// doc/drawing_TrianglePETRing.png
/*
1.
In ANY triangle with sides with lengths a, b, c, and an angle alpha between b and c (opposite a), the following is true:
    a^2 = b^2 + c^2 - 2 * b * c * cos(alpha)

2.
In our problem we have:
   - A source (distance S from source to origin: vec_o - vec_s)
   - A hit1 (distance D1 from hit-1 to source: vec_s - vec_hit_1)
   - A radius R
   - The angle alpha (this is the tricky bit): cos(alpha) = (vec_s - vec_hit_1) * (vec_o - vec_s)/(lengths of both vectors)
     We do this, with the call:
        double cosAngle = CalculateCosAngleBD( hit1ToSrc, srcToO );

3.
   - What we want to know is F: the distance from source to hit-2

4.
   - Triangle formula:
        R^2 = F^2 + S^2 - 2 * F * S * cos(alpha)
   - Rewrite:
        F^2 - 2 * S * cos(alpha) * F + (S^2 - R^2) = 0

5.
   - General Quadratic formula:
        If: a*x^2 + b*x + c = 0
        Then: x = [-b +- SQRT(b^2 - 4ac)] / 2a
   - In our case:
        a = 1
        b = - 2 * S * cos(alpha)
        c = (S^2 - R^2)
   - Call solveQuadraticEquation gives solution_1 and solution_2 (should be larger than 0 since are talking about a distance)

6.
   - The location of hit-2:
        hit-2 = vec_s + (vec_s - vec_hit_1)/length_of_vector * F

*/

// triangle with sides A, B, C and angle alpha_opposite_a (between sides B and C, opposite side A).
// unknown: length B
//
int CalculateLengthBOfTriangle(const double& in_A, const double& in_C, const double& in_cosAngle_A,
                           double& out_solution1, double& out_solution2 )
{
    // coefficients of cuadratic formula: a * x^2 + b*x + c = 0;
    double a = 1.0;
    double b = -2.0 * in_C * in_cosAngle_A;
    double c = in_C*in_C - in_A*in_A;
    QuadraticFactors qfactors(a, b, c);

    //  cout << "quad f, a, b, c: " << a << ", " << b << ", " << c << endl;
    int numsols = solveQuadraticEquation(qfactors, out_solution1, out_solution2);
    //  cout << "numsols: " << numsols << ": " << out_solution1 << " " << out_solution2 << endl;

    return numsols;
}

