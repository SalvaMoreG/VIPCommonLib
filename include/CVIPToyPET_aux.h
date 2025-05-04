
#pragma once

#ifndef _CVIPToyPET_aux___
#define _CVIPToyPET_aux___

#include "CVIPUtils.h"
#include "CVIP3Vector.h"

// triangle with sides A, B, C and angle alpha_opposite_a (between sides B and C, opposite side A).
// unknown: length B
//
int CalculateLengthBOfTriangle(const double& in_A, const double& in_C, const double& in_cosAngle_A,
                           double& out_solution1, double& out_solution2 );
/*
In ANY triangle with sides with lengths A, B, C, and an angle alpha_opposite_A between B and C (opposite A),
the following is true:
    A^2 = B^2 + C^2 - 2 * B * C * cos(alpha_opposite_A)
Hence:
    B^2 - 2 * C * cos(alpha_opposite_A) * B + (C^2 - A^2) = 0
Quadratic equation:
    unknown B -> x
        a = 1
        b = - 2 * C * cos(alpha_opposite_A)
        c = (C^2 - A^2)
*/


// returns true if a new hit position was found after scattering
bool CalculateScatteredTrack( const double& in_radiusPET, const double& in_radiusSphere
           , const C3Vector& in_sourcePosition, C3Vector& io_hitPosition, C3Vector& out_scatterPosition
           , double& io_newEnergy );

bool CalculateSpatialSmearedHitPosition( const double& in_radiusPET, C3Vector& io_hitPosition );

void AuxFinalizeH();

#endif
