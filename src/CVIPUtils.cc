

#include "VIPconstants.h"
#include "CVIPUtils.h"
#include "CVIP3Vector.h"

#include "CVIPRandom.h"

#include <cassert>
#include <sys/stat.h>

using namespace std;

// ****************************************************************************

bool doubleEquals(double left, double right, double epsilon)
{
	// return (fabs(left - right) <= epsilon * fabs(left));
	/*
	// TODO: check this a whole lot more!!!
	//  The problem is that in the case of matrix multiplication, we can have:
	// 		left = x
	//		right = 0
	// 	and so, the comparison " x <= epsilon * x" is always false....
	*/
	return (fabs(left - right) <= epsilon);
}

// ****************************************************************************************************

bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

/*
// TODO: which version is better????
bool fileExists(const std::string& in_filename)
{
    ifstream infile(in_filename.c_str());
    return infile.good();
}
*/

// ****************************************************************************************************

void
PrintAxis( AXIS in_axis )
{
	if (in_axis == Axis_X)
		cout << " Axis_X ";
	else if (in_axis == Axis_Y)
		cout << " Axis_Y ";
	else if (in_axis == Axis_Z)
		cout << " Axis_Z ";
	else
		assert( false );
}

// ****************************************************************************

void
WriteString(std::ostream& io_os, const std::string& in_string)
{
	int size = in_string.length();
	io_os.write( (char*) &size, sizeof(size) );
	io_os.write( &in_string[0], size );
}

// ****************************************************************************

void
ReadString(std::istream& in_is, std::string& io_string)
{
   size_t length;
   in_is.read((char*)&length, sizeof(size_t));
   in_is.width(length);
   in_is >> io_string;
}

// ****************************************************************************

void bubbleSort( int size, double* x )
{
   int n = size;
   while (n != 0)
   {
      int newn = 0;
      for (int i = 1; i < n; i++)
      {
         if (x[i-1] > x[i])
         {
            double tmp = x[i];
            x[i] = x[i-1];
            x[i-1] = tmp;
            newn = i;
         }
      }
      n = newn;
   }
}

// ****************************************************************************

// TODO: there are other ways to implement this (more exact to double rounding errors...)?
int
solveQuadraticEquation(const QuadraticFactors& in_quadraticFactors, double& out_x1, double& out_x2)
{
	out_x1 = 0;
	out_x2 = 0;

	double discriminant = (in_quadraticFactors.GetbFactor() * in_quadraticFactors.GetbFactor())
				- (4 * in_quadraticFactors.GetaFactor() * in_quadraticFactors.GetcFactor());

	if (    in_quadraticFactors.GetaFactor() == 0 && in_quadraticFactors.GetbFactor() == 0
		 && in_quadraticFactors.GetcFactor() == 0)
	{
		return 0;
	}

	// If a == 0, equation becomes linear (cannot divide by 0)
	if ( in_quadraticFactors.GetaFactor() == 0 )
	{
		if (in_quadraticFactors.GetbFactor() != 0)
		{
			out_x1 = -1.0 * in_quadraticFactors.GetcFactor()/in_quadraticFactors.GetbFactor();
			return 1;	// 1 solution
		}
		else
			return 0;	// if b = 0, there are no solutions (infinite solutions)
	}

	// old stuff
    int numSolutions = 0;
    if (discriminant >= 0)
    {
        out_x1 = ((-1 * in_quadraticFactors.GetbFactor()) + (sqrt(discriminant)))/(2 * in_quadraticFactors.GetaFactor());
        numSolutions++;

        if (discriminant > 0)
        {
            out_x2 = ((-1 * in_quadraticFactors.GetbFactor()) - (sqrt(discriminant)))/(2 * in_quadraticFactors.GetaFactor());
            numSolutions++;
        }
	}
    return numSolutions;
}

// ****************************************************************************

double GetPhiFromXandY(const double& in_X, const double& in_Y)
{
	double out_phi;
	if (in_X == 0)
	{
		if (in_Y > 0) out_phi = 0.5 * kPI;
		else out_phi = 1.5 * kPI;
	}
	else if (in_Y == 0)
	{
		if (in_X > 0) out_phi = 0.0;
		else out_phi = kPI;
	}
	else
	{
		double quotient = fabs(in_Y / in_X);
		if (in_X > 0 && in_Y > 0)
		{
			out_phi = atan(quotient);
		}
		else if (in_X < 0 && in_Y < 0)
		{
			out_phi = kPI + atan(quotient);
		}
		else if (in_X < 0 && in_Y > 0)
		{
			out_phi = kPI - atan(quotient);
		}
		else if (in_X > 0 && in_Y < 0)
		{
			out_phi = 2 * kPI - atan(quotient);
		}
	}

	return out_phi;
}

// ****************************************************************************

bool
RayBoxIntersection(const C3Vector& in_origin, const C3Vector& in_direction
				 , const C3Vector& lowerEdge, const C3Vector& upperEdge
				 , double& io_factor_enter, double& io_factor_exit )
{
/*
	Ray/box intersection using the Smits' algorithm

	Input:
		origin.
		direction.
		box = (vmin,vmax)
	Output:
		flag: (0) Reject, (1) Intersect.
		tmin: distance from the ray origin.
	Author:
		Jesus Mena
*/

/*
// ALSO:
	* Ray-box intersection using IEEE numerical properties to ensure that the
	* test is both robust and efficient, as described in:
	*
	*      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
	*      "An Efficient and Robust Ray-Box Intersection Algorithm"
	*      Journal of graphics tools, 10(1):49-54, 2005
	*
	// NOTE 1: Division by 0 (or -0) could be taken care of by IEEE standards "inf" and "-inf"
	// 		because in principle C++ can correctly compare (-)inf with valid numbers
	// 		however, I am not sure this is compiler- and/or system independent
	// NOTE 2: I took for the minimum value for tmin = 0 and the maximum value for tmax = 1.0
	// 		This is because the point "t" should always lie in between point 1 and point 2, never outside them.
*/


	double epsilon = 0.0001;
	bool xValid = ( doubleEquals(in_direction.GetX(), 0.0, epsilon) ) ? false : true;
	bool yValid = ( doubleEquals(in_direction.GetY(), 0.0, epsilon) ) ? false : true;
	bool zValid = ( doubleEquals(in_direction.GetZ(), 0.0, epsilon) ) ? false : true;

/*
cout << "xyz valids: " << xValid << " " << yValid << " " << zValid << endl;
cout << "in_origin: " << in_origin << endl;
*/

	if (!xValid && !yValid && !zValid)
		return false;

	if (!xValid && ( in_origin.GetX() < lowerEdge.GetX() || in_origin.GetX() > upperEdge.GetX() ))
    {
		return false;
    }
	if (!yValid && ( in_origin.GetY() < lowerEdge.GetY() || in_origin.GetY() > upperEdge.GetY() ))
    {   
		return false;
    }
	if (!zValid && ( in_origin.GetZ() < lowerEdge.GetZ() || in_origin.GetZ() > upperEdge.GetZ() ))
    {
		return false;
    }

// cout << "I am still here" << endl;

	double tmin = 0.0, tmax = 1.0;
	if (xValid)
	{
		if (in_direction.GetX() > 0)
		{
			tmin = (lowerEdge.GetX() - in_origin.GetX()) / in_direction.GetX();
            /*
            cout <<"tmin calc(1), lowEdgeX: "<< lowerEdge.GetX()<<" origin.X: "<< in_origin.GetX() 
                 << " d.X: " << in_direction.GetX() 
                << " lX - oX: " << (lowerEdge.GetX() - in_origin.GetX())
                << " tmin: " << tmin << endl;
                */
			tmax = (upperEdge.GetX() - in_origin.GetX()) / in_direction.GetX();
            /*
            cout <<"tmax calc(1), uppEdgeX: "<< upperEdge.GetX()<<" origin.X: "<< in_origin.GetX() 
                 << " d.X: " << in_direction.GetX() 
                << " uX - oX: " << (upperEdge.GetX() - in_origin.GetX())
                << " tmax: " << tmax << endl;
                */
		}
		else if (in_direction.GetX() < 0)
		{
			tmin = (upperEdge.GetX() - in_origin.GetX()) / in_direction.GetX();
            /*
            cout << "tmin calc(2), uppEdgeX: " << upperEdge.GetX() << " origin.X: " << in_origin.GetX() 
                 << " d.X: " << in_direction.GetX() 
                << " uX - oX: " << (upperEdge.GetX() - in_origin.GetX())
                << " tmin: " << tmin << endl;
                */
			tmax = (lowerEdge.GetX() - in_origin.GetX()) / in_direction.GetX();
            /*
            cout << "tmax calc(2), lowEdgeX: " << lowerEdge.GetX() << " origin.X: " << in_origin.GetX() 
                 << " d.X: " << in_direction.GetX() 
                << " lX - oX: " << (lowerEdge.GetX() - in_origin.GetX())
                << " tmax: " << tmax << endl;
                */
		}
		//    cout << "X: " << in_direction.GetX() << " tmin: " << tmin << " tmax: " << tmax << endl;
	}

	double tymin = 0.0, tymax = 1.0;
	if (yValid)
	{
		if (in_direction.GetY() > 0)
		{
			tymin = (lowerEdge.GetY() - in_origin.GetY()) / in_direction.GetY();
			tymax = (upperEdge.GetY() - in_origin.GetY()) / in_direction.GetY();
		}
		else if (in_direction.GetY() < 0)
		{
			tymin = (upperEdge.GetY() - in_origin.GetY()) / in_direction.GetY();
			tymax = (lowerEdge.GetY() - in_origin.GetY()) / in_direction.GetY();
		}
		//    cout << "Y: " << in_direction.GetY() << " tymin: " << tymin << " tymax: " << tymax << endl;
	}

    if ( (tmin > tymax) || (tymin > tmax) )
	{
        // Nothing found...
		// cout << "Y ts: " << tmin << " " << tmax << " tys: " << tymin << " " << tymax << endl;
        tmin = -1;
    	return false;
	}

    if (tymin > tmin)
        tmin = tymin;

	if (tymax < tmax)
        tmax = tymax;

	double tzmin = 0.0, tzmax = 1.0;
	if (zValid)
	{
		if (in_direction.GetZ() >= 0)
		{
			tzmin = (lowerEdge.GetZ() - in_origin.GetZ()) / in_direction.GetZ();
			tzmax = (upperEdge.GetZ() - in_origin.GetZ()) / in_direction.GetZ();
		}
		else
		{
			tzmin = (upperEdge.GetZ() - in_origin.GetZ()) / in_direction.GetZ();
			tzmax = (lowerEdge.GetZ() - in_origin.GetZ()) / in_direction.GetZ();
		}
		//    cout << "Z: " << in_direction.GetZ() << " tzmin: " << tzmin << " tzmax: " << tzmax << endl;
	}

    if ((tmin > tzmax) || (tzmin > tmax))
	{
        // Nothing found...
        // cout << "Z ts: " << tmin << " " << tmax << " tzs: " << tzmin << " " << tzmax << endl;
        tmin = -1;
    	return false;
	}

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

	io_factor_enter = (tmin < tmax) ? tmin : tmax;
	io_factor_exit  = (tmin < tmax) ? tmax : tmin;
/*
cout << "final RayBoxfactors tmin: " << tmin << "  " << tmax << endl;
cout << "final RayBoxfactors: " << io_factor_enter << "  " << io_factor_exit << endl;
*/
	// Make sure the factors are in between 0.0 and 1.0 (if not, the ray does not go through the FOV)

// cout << "fac_enter: " << io_factor_enter <<  " factor-exit: " << io_factor_exit << endl;

	// return ( io_factor_enter < 1.0 && io_factor_exit > 0.0 ) && ( io_factor_enter > 0.0 && io_factor_exit < 1.0 );
    
	return ( io_factor_enter > 0.0 && io_factor_enter < 1.0 
	       && io_factor_exit > 0.0 && io_factor_exit < 1.0 );
}

// ****************************************************************************

/*
// =================================================================================
 double NextRayIntersection()

 Input:
	- origin of the "ray"
	- direction of the "ray"
	- 3-vector with lower bounds of current FOV bin
	- 3-vector with upper bounds of current FOV bin
 Output:
	- fraction of direction to get to new FOV bin

 Method:
	Looks for the next bounds (lower or upper) in the FOV in the X, Y or Z direction,
	moving along the "direction" of the ray. The smallest fraction of the "direction"
	starting from the current position gives the position of the next intersection.
*/

double
NextRayIntersection( const C3Vector& in_origin, const C3Vector& in_direction
	, const C3Vector& in_minBounds, const C3Vector& in_maxBounds )
{
	double factor = 1.0;

	// newbounds = current + factor * direction
	// factor * direction = newbounds - current
	// factor = ( newbounds - current ) / direction
	if (in_direction.GetX() > 0)
		factor = (in_maxBounds.GetX() - in_origin.GetX())/in_direction.GetX();
	else if (in_direction.GetX() < 0)
		factor = (in_minBounds.GetX() - in_origin.GetX())/in_direction.GetX();

	// same for y
	double factorB = 1;
	if (in_direction.GetY() > 0)
		factorB = (in_maxBounds.GetY() - in_origin.GetY())/in_direction.GetY();
	else if (in_direction.GetY() < 0)
		factorB = (in_minBounds.GetY() - in_origin.GetY())/in_direction.GetY();

	//
	factor = (factor < factorB) ? factor : factorB;

	// same for z
	if (in_direction.GetZ() > 0)
		factorB = (in_maxBounds.GetZ() - in_origin.GetZ())/in_direction.GetZ();
	else if (in_direction.GetZ() < 0)
		factorB = (in_minBounds.GetZ() - in_origin.GetZ())/in_direction.GetZ();

	//
	factor = (factor < factorB) ? factor : factorB;

	// TODO (TOTHINKABOUT:)
	// HMM, as a matter of fact, this means (ix++, iy, iz) OR (ix, iy++, iz) OR (ix, iy, iz++)

	return factor;

}

// ****************************************************************************

bool IsInsideBounds(const C3Vector& in_minBounds, const C3Vector& in_maxBounds, const C3Vector& in_position)
{
	bool insideBounds = (  in_position.GetX() >= in_minBounds.GetX()
						&& in_position.GetX() <= in_maxBounds.GetX()
						&& in_position.GetY() >= in_minBounds.GetY()
						&& in_position.GetY() <= in_maxBounds.GetY()
						&& in_position.GetZ() >= in_minBounds.GetZ()
						&& in_position.GetZ() <= in_maxBounds.GetZ() );
	return insideBounds;
}

// ****************************************************************************

void
VIPUtils::GetRandomSolidAngle(double& io_phi, double& io_theta, bool doFourPi)
{
	CVIPRandomUniform uniformDist(0.0, 1.0);
	io_phi = 2.0 * kPI * uniformDist.GetNewValue();

	double cosTheta = 1.1;
	while ( fabs(cosTheta) > 1.0 )
	{
		cosTheta = uniformDist.GetNewValue();
            // generate flat in cos(theta), with cos(th) = [0, 1]
			// (If we only care for events that get into detector (i.e., angle < 90 degrees))
		// Else, if generating over 4 pi:
		if (doFourPi)
		{
			cosTheta = 1 - 2*cosTheta;				// cos = [-1, 1]
		}
	}
	io_theta = acos( cosTheta );  		// acos is inverse cosine
}

// ****************************************************************************

void
VIPUtils::GetRandomPositionOnSolidAngle( const double& in_phi, const double& in_theta
				, const double& in_zMin, const double& in_zSize, const C3Vector& in_originPosition
				, C3Vector& out_randomPosition )
{
	double px = sin(in_theta) * cos(in_phi);
	double py = sin(in_theta) * sin(in_phi);
	double pz = cos(in_theta);

	CVIPRandomUniform uniformDist(0, in_zSize);
	double posZ = in_zMin + uniformDist.GetNewValue();

	double factor = (posZ - in_originPosition.GetZ()) / pz;

	double posX = in_originPosition.GetX() + factor * px;
	double posY = in_originPosition.GetY() + factor * py;

	out_randomPosition.Set(posX, posY, posZ);
}

// ****************************************************************************

// Nov 26
// By using different spheres (with R = distance from origin to sphere surface position), the "ray/area" density can be easily calculated!!
// Using the simple fact that the density of points on a surface is: D = N_generated / radius^2
//
double
VIPUtils::GetPseudoSolidAngleProbability( const C3Vector& in_originCoords, const C3Vector& in_surfaceCoords )
{
	// 1. Get vector from source (or fov-voxel-centre) to hit-position in scatterer
	C3Vector radialVector = in_surfaceCoords - in_originCoords;
	double zDistance = radialVector.GetZ();

	// density 1, which is the maximum density
	// density_1 = N / area_1
	double area_1 = zDistance * zDistance; 		// area = 4 pi r^2

	// angle
	double transverse = sqrt( radialVector.GetX()*radialVector.GetX() + radialVector.GetY()*radialVector.GetY() );
	double tanTh = transverse / zDistance;
	double theta = atan(tanTh);

	// density 2 = N / area_2
	double radius_2 = zDistance / cos( theta );
	double area_2 = radius_2 * radius_2;

	// probability = density_2 / density_1 = area_1 / area_2; (so that fraction = [0, 1])
	double fraction = area_1 / area_2;
	return fraction;
}

// ****************************************************************************

double
VIPUtils::Gauss( const C3Vector& in_mean, const double& in_sigma, const C3Vector& in_coords, int in_dim )
{
	double factor = sqrt(2.0 * kPI) * in_sigma;
	double term1 = factor;
	for (int i = 0; i < in_dim; i++)
		term1 = term1*factor;

	double sigma2 = in_sigma*in_sigma;
	C3Vector difvec = in_coords-in_mean;
	double term2 = -1.0 * (difvec*difvec)/(2.0 * sigma2);
	double gauss = term1 * exp(term2);

	return gauss;
}


// ===========================================================

double
VIPUtils::minAbsDifferenceTwoAnglesRad( const double& in_angle1, const double& in_angle2 )
{
    double diff = abs(in_angle1 - in_angle2); 
    if (diff > kPI)
        diff = 2.0*kPI - diff;
    return diff;    
}

// ****************************************************************************

bool
VIPUtils::GetGeomAngleDegrees(const C3Vector& v1, const C3Vector& v2, double& out_angleGeomDegrees)
{
   out_angleGeomDegrees = -1000;

   double vectorproduct = v1*v2;
   double len1 = v1.GetLength();
   double len2 = v2.GetLength();
   if (len1 > 0 && len2 > 0)
   {
      double cosTheta = vectorproduct/(len1*len2);
      if (cosTheta > -1.0 && cosTheta < 1.0 )
      {
	     out_angleGeomDegrees = acos( cosTheta ) * 180.0/kPI;
	     // cout << "geometric angle: " << out_thetaGeom << endl;
	     return true;
      }
   }
   return false;
}

// ======================================================================

bool
VIPUtils::GetPositionAngle2D( const double& in_x, const double& in_y, double& out_angleGeom)
{
	bool debug(false);
	if (in_x == 0)
	{
		out_angleGeom = (in_y > 0) ? 0.5 * kPI : 1.5 * kPI;
	}
	else if (in_y == 0)
	{
		out_angleGeom = (in_x > 0) ? 0.0 : kPI;
	}
	else if (in_x > 0 && in_y > 0)
	{
		out_angleGeom = atan(in_y/in_x);
		if (debug) cout << "(" << in_x << ", " << in_y << "): "
			 << in_y << " / " << in_x << " : " << in_y/in_x
			 << " angle: " << out_angleGeom << " = " << out_angleGeom*180/kPI << endl;
	}
	else if (in_x < 0 && in_y > 0)
	{
		out_angleGeom = 0.5*kPI + atan(fabs(in_x/in_y));
		if (debug) cout << "(" << in_x << ", " << in_y << "): "
			 << " fabs(x/y): " << fabs(in_x/in_y) << " atan: " << atan(fabs(in_x/in_y)) << " + 0.5PI "
			 << " angle: " << out_angleGeom << " = " << out_angleGeom*180/kPI << endl;

	}
	else if (in_x < 0 && in_y < 0)
	{
		out_angleGeom = kPI + atan(fabs(in_y/in_x));
		if (debug) cout << "(" << in_x << ", " << in_y << "): "
			 << " fabs(y/x): " << fabs(in_y/in_x) << " atan: " << atan(fabs(in_y/in_x)) << " + PI "
			 << " angle: " << out_angleGeom << " = " << out_angleGeom*180/kPI << endl;
	}
	else if (in_x > 0 && in_y < 0)
	{
		out_angleGeom = 2*kPI - atan(fabs(in_y/in_x));
		if (debug) cout << "(" << in_x << ", " << in_y << "): "
			 << " fabs(y/x): " << fabs(in_y/in_x) << " atan: " << atan(fabs(in_y/in_x)) << " 2PI-atan "
			 << " angle: " << out_angleGeom << " = " << out_angleGeom*180/kPI << endl;
	}
	else
	{
		out_angleGeom = 0.0;
		cout << "VIPUtils::GetPositionAngle2D, ERROR: " << in_x << " " << in_y << endl;
		return false;
	}
	return true;
}


// ======================================================================

bool
VIPUtils::GetComptonAngleDegrees(const double& E1_keV, const double& Etot_keV, double& out_angleComptonDegrees )
{
   out_angleComptonDegrees = -1000;

   if (Etot_keV > 0 && (Etot_keV - E1_keV) > 0)
   {
       double coscompton = 1.0 - mass_electron_keV*(1.0/(Etot_keV - E1_keV) - 1.0/Etot_keV);
       if (fabs(coscompton) <= 1)
       {
          out_angleComptonDegrees = acos(coscompton) * 180.0/kPI;
          return true;
       }
   }
   return false;
}

bool
VIPUtils::GetComptonEnergyFromAngleRad(const double& in_angleComptonRadians, const double& Etot_keV, double& out_E1_keV )
{
	double cosTh = cos(in_angleComptonRadians);
	double tmp = 1.0 + mass_electron_keV / (Etot_keV * (1.0 - cosTh)) ;

	out_E1_keV = Etot_keV / tmp;

	// Since all angles are permitted (0 - 180 degrees; higher than 180 degrees gives same cosine as smaller than 180 degrees)
	// this function always returns "true"
	return true;
}

// ======================================================================

double
VIPUtils::GetDistance(const C3Vector& v1, const C3Vector& v2)
{
	C3Vector tmp = v2 - v1;
	return tmp.GetLength();
}

// ======================================================================

void
VIPUtils::GetFileNames ( const string& in_string, vector<string>& out_fileNames ) // COMMA separated substrings
{
	/* // No spaces possible in string"
	string tmpStr( in_string );
 	std::string::iterator end_pos = std::remove(tmpStr.begin(), tmpStr.end(), ' ');
	tmpStr.erase(end_pos, tmpStr.end());
	*/

	stringstream ss( in_string );
	out_fileNames.clear();
	while( ss.good() )
	{
   		string substr;
   		getline( ss, substr, ',' );
   		out_fileNames.push_back( substr );
	}
}


// ===================================================================

double
VIPUtils::PolynomialGrade3( double* in_x, double* in_par )
{
	double x = in_x[0]; // the fit variable is given as a pointer in ROOT
	double results = in_par[0] + in_par[1]*x + in_par[2]*x*x + in_par[3]*x*x*x;
	return results;
}

// ===================================================================

double
VIPUtils::PolynomialGrade4( double* in_x, double* in_par )
{
	double x = in_x[0]; // the fit variable is given as a pointer in ROOT
	double results = in_par[0] + in_par[1]*x + in_par[2]*x*x + in_par[3]*x*x*x + in_par[4]*x*x*x*x;
	return results;
}

// ===================================================================

// Explanation:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
// 
bool
VIPUtils::CalculateIntersectionDistance(const C3Vector& in_X1, const C3Vector& in_X2
        , const C3Vector in_Centre, const double& in_radius, double& out_distance)
{
    C3Vector LOR = in_X2 - in_X1;
    C3Vector unitLOR(LOR);
    unitLOR.Normalize();
    
    C3Vector XtoC = in_Centre - in_X1;
    
    double angleLOR_XtoC = XtoC.GetScalarProductAngleRadians(LOR);
    double t_ca = (angleLOR_XtoC < 0.0001) ? XtoC.GetLength() : XtoC * unitLOR;
    
    bool debug(false);
    bool failed(false);
    if (t_ca < 0) 
    {
        if (debug)
        {
            cout << "solution not found" << endl;
            cout << "X1 " << in_X1 << " X2: " << in_X2 
                << " LOR: " << LOR << " X1toC: " << XtoC << endl;
            cout << "t_ca = " << t_ca << endl;
        }
        failed = true;
    }
    
    double XtoC_len, CtoLOR_len; 
    if (!failed)
    {
        XtoC_len = XtoC.GetLength();
        CtoLOR_len = sqrt(XtoC_len * XtoC_len - t_ca * t_ca);
        if (CtoLOR_len < 0)
        {
            if (debug)
            {
                cout << "solution not found (2)" << endl;
                cout << "X1 " << in_X1 << " X2: " << in_X2 
                     << " LOR: " << LOR << " X1toC: " << XtoC << endl;
                cout << "CtoLOR_len = " << CtoLOR_len << endl;
            }
            failed = true;
        }
    }
    
    if (!failed)
    {
        if (CtoLOR_len > in_radius)
        {
            if (debug)
            {
                cout << "LOR does not touch sphere" << endl;
                cout << "X1 " << in_X1 << " X2: " << in_X2 
                    << " LOR: " << LOR << " X1toC: " << XtoC << endl;
                cout << "CtoLOR_len = " << CtoLOR_len << endl;
            }
            failed = true;
        }
    }

    if (!failed)
    {
        double t_hc = sqrt(in_radius*in_radius - CtoLOR_len*CtoLOR_len);
        double t0 = t_ca - t_hc;
        double t1 = t_ca + t_hc;
        
        C3Vector P0 = in_X1 + unitLOR * t0;
        C3Vector P1 = in_X1 + unitLOR * t1;
        out_distance = (P1 - P0).GetLength();
    }

    return (!failed);
}
