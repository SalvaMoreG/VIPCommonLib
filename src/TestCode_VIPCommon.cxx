 
#include <iostream>
#include <fstream>

#include "CVIPGaussianElimination.h"
#include "VIPconstants.h"
#include "CVIPUtils.h"
#include "CVIP3Vector.h"
#include "CVIPFieldOfView.h"
#include "CVIPWildermanUtils.h"

#include "CVIPImageKernel.h"
#include "CVIPGeometry.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1D.h"

#include <vector>
#include <cassert>

#include "CVIPRandom.h"

using namespace std;

void TestWildermanMarch();
void TestMarchCircle();
void TestFromConeEquationToQuadraticFactor();

void TestRotateFieldOfView();

void TestImageKernel();

void TestVectorsAndMatrices();

void TestVipGeometry();

void TestComptonEdge();

void TestSolveSetOfEquations_grade3();
void TestSolveSetOfEquations_grade4();
void CalculateP0_P1_Approx(const double* in_x, const double* in_y, double& out_p0, double& out_p1);
// -------------

void CreateSimpleFOV( CVIPFieldOfView& io_fieldOfView );
void CreateSimpleCircle( const double& in_sliceZ, const double& in_circleradius
	, const double& origin_x, const double& origin_y, C3Vector& io_origin
	, const double& in_axis_x, const double& in_axis_y, C3Vector& io_coneAxis, double& io_comptonAngle);

void TestVectorSubtraction();

void TestVIPRandom();
void TestVIPRandom2();
void AuxGetTheRandoms( double& aUniform, double& aGauss);

// ***********************************************************************************************

int main()
{
	// TestWildermanMarch();

	// TestMarchCircle();

	// TestFromConeEquationToQuadraticFactor();

	// TestRotateFieldOfView();

	// TestImageKernel();

	// TestVectorsAndMatrices();

	// TestVipGeometry();
	
	// TestComptonEdge();

	//	TestSolveSetOfEquations_grade3();
	//	TestSolveSetOfEquations_grade4();

	//  TestVectorSubtraction();

    TestVIPRandom();
	TestVIPRandom2();
}

// ***********************************************************************************************


void
TestWildermanMarch()
{
	cout << endl;
	cout << "*** TestWildermanMarch() ***" << endl;

	CVIPFieldOfView fieldOfView;
	CreateSimpleFOV( fieldOfView );

	cout << endl;
	cout << "Circle, center at (-5, 5), radius 3.5 " << endl;
	{
		// easy circle (centre at -5, 5)
		double circleradius = 3.5;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, -5, 5, origin, 0, 0, coneAxis, comptonAngle );

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Circle, center at (-5, 5), radius 1 " << endl;
	{
		// easy circle (centre at -5, 5), radius 1
		double circleradius = 1.0;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, -5, 5, origin, 0, 0, coneAxis, comptonAngle );

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Circle, center at (-1.5, 1.5), radius 3.5 " << endl;
	{
		// easy circle  (centre at -1.5, 1.5)
		double circleradius = 3.5;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;

		CreateSimpleCircle( sliceZ, circleradius, -1.5, 1.5, origin, 0, 0, coneAxis, comptonAngle );

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Circle, center at (-1.5, 1.5), radius 2 -> scenario 1 " << endl;
	{
		// easy circle (centre at -1.5, 1.5), radius 2
		//  --> scenario 1
		double circleradius = 2.0;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, -1.5, 1.5, origin, 0, 0, coneAxis, comptonAngle );

		// Get the ellips parameters
		EllipsParameters ellipsParameters;
		CalculateEllipsFactors( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, ellipsParameters );
		cout << "ellips parameters, A: " << ellipsParameters.A << " B: " << ellipsParameters.B << " C: " << ellipsParameters.C
			 << " D: " << ellipsParameters.D << " E: " << ellipsParameters.E << " F: " << ellipsParameters.F << endl;
		cout << "more ellips parameters, axis_a: " << ellipsParameters.axis_a << " axis_b: " << ellipsParameters.axis_b 
			   << " centre_x: " << ellipsParameters.centre_x << " centre_y: " << ellipsParameters.centre_y << endl;

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Circle, center at ( 0, 0), radius 3.2 -> scenario 1 " << endl;
	{
		// easy circle (centre at 0, 0), radius 2
		//  --> scenario 1
		double circleradius = 3.2;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, 0, 0, origin, 0, 0, coneAxis, comptonAngle );

		// Get the ellips parameters
		EllipsParameters ellipsParameters;
		CalculateEllipsFactors( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, ellipsParameters );
		cout << "ellips parameters, A: " << ellipsParameters.A << " B: " << ellipsParameters.B << " C: " << ellipsParameters.C
			 << " D: " << ellipsParameters.D << " E: " << ellipsParameters.E << " F: " << ellipsParameters.F << endl;
		cout << "more ellips parameters, axis_a: " << ellipsParameters.axis_a << " axis_b: " << ellipsParameters.axis_b 
			   << " centre_x: " << ellipsParameters.centre_x << " centre_y: " << ellipsParameters.centre_y << endl;

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Circle, center at ( 0, 0), radius 2 -> scenario 1 " << endl;
	{
		// easy circle (centre at 0, 0), radius 2.
		//  --> scenario 1
		double circleradius = 2.;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, 0, 0, origin, 0, 0, coneAxis, comptonAngle );

		// Get the ellips parameters
		EllipsParameters ellipsParameters;
		CalculateEllipsFactors( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, ellipsParameters );
		cout << "ellips parameters, A: " << ellipsParameters.A << " B: " << ellipsParameters.B << " C: " << ellipsParameters.C
			 << " D: " << ellipsParameters.D << " E: " << ellipsParameters.E << " F: " << ellipsParameters.F << endl;
		cout << "more ellips parameters, axis_a: " << ellipsParameters.axis_a << " axis_b: " << ellipsParameters.axis_b 
			   << " centre_x: " << ellipsParameters.centre_x << " centre_y: " << ellipsParameters.centre_y << endl;

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}

	cout << endl;
	cout << "Ellips, derived from center at (0, 0), axis direction = (x, y, z), radius = 2 -> scenario 1 " << endl;
	// See: http://www.maa.org/joma/volume8/kalman/QuadForm.html
	{
		// easy circle (centre at -2, 0), radius 2.
		//  --> scenario 1
		double circleradius = 2.;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, -2, 0, origin, 3, 0, coneAxis, comptonAngle );

		// Get the ellips parameters
		EllipsParameters ellipsParameters;
		CalculateEllipsFactors( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, ellipsParameters );
		cout << "ellips parameters, A: " << ellipsParameters.A << " B: " << ellipsParameters.B << " C: " << ellipsParameters.C
			 << " D: " << ellipsParameters.D << " E: " << ellipsParameters.E << " F: " << ellipsParameters.F << endl;
		cout << "more ellips parameters, axis_a: " << ellipsParameters.axis_a << " axis_b: " << ellipsParameters.axis_b 
			   << " centre_x: " << ellipsParameters.centre_x << " centre_y: " << ellipsParameters.centre_y << endl;

		// Do the Wilderman March
		std::vector<unsigned int> fovidxVec;
		int scenario = WildermanMarch( sliceZ, origin, coneAxis, comptonAngle, fieldOfView, fovidxVec );
		cout << "SCENARIO: " << scenario << endl;

		ShowFovVector( fieldOfView, fovidxVec );
	}
}

// -----------------------------------------------------------------------------------------

void
TestMarchCircle()
{
	cout << endl;
	cout << "*** TestMarchCircle() ***" << endl;

	CVIPFieldOfView fieldOfView;
	CreateSimpleFOV( fieldOfView );

	{
		// easy circle
		double circleradius = 3.5;
		double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
		C3Vector origin, coneAxis;
		double comptonAngle, dummy = 0.0;
		CreateSimpleCircle( sliceZ, circleradius, -5, 5, origin, 0, 0, coneAxis, comptonAngle );
	
		// Find solutions....
		// find an intersection on the X axis
		C3Vector edge( fieldOfView.GetLowerEdge().GetX(), dummy, sliceZ );
		C3Vector sol1, sol2;
		QuadraticFactors quadraticFactors;
		int num = FindIntersectionsInFOV( origin, coneAxis, comptonAngle, edge, Axis_Y
						, fieldOfView.GetLowerEdge(), fieldOfView.GetUpperEdge(), sol1, sol2, quadraticFactors );
		assert( num == 1);
		C3Vector solxmin( sol1 );
	
		// find an intersection on the Y axis
		edge.Set( dummy, fieldOfView.GetUpperEdge().GetY(), sliceZ );
		num = FindIntersectionsInFOV( origin, coneAxis, comptonAngle, edge, Axis_X
						, fieldOfView.GetLowerEdge(), fieldOfView.GetUpperEdge(), sol1, sol2, quadraticFactors );
		assert( num == 1);
		C3Vector solymax( sol1 );
	
		cout << endl;
	
		// march function
		std::vector<unsigned int> fovidxVec;
		// From xmin to ymax
		SEARCH_DIRECTION searchDirection( SEARCH_NEXT_X );
		March( origin, coneAxis, comptonAngle, solxmin, solymax, fieldOfView, fovidxVec, searchDirection );
	
		ShowFovVector( fieldOfView, fovidxVec );
	
		// From ymax to xmin
		searchDirection = SEARCH_PREV_Y ;
		March( origin, coneAxis, comptonAngle, solymax, solxmin, fieldOfView, fovidxVec, searchDirection );
	
		ShowFovVector( fieldOfView, fovidxVec );
	}
}

// -----------------------------------------------------------------------------------------

void
TestFromConeEquationToQuadraticFactor()
{
	cout << endl;
	cout << "*** TestFromConeEquationToQuadraticFactor() ***" << endl;

	CVIPFieldOfView fieldOfView;
	CreateSimpleFOV( fieldOfView );

	// easy circle
	double circleradius = 3.5;
	double sliceZ = 0.5 * (fieldOfView.GetUpperEdge().GetZ() + fieldOfView.GetLowerEdge().GetZ());
	C3Vector origin, coneAxis;
	double comptonAngle;

	CreateSimpleCircle( sliceZ, circleradius, -5, 5, origin, 0, 0, coneAxis, comptonAngle );

	double dummy;
	QuadraticFactors quadraticFactors, quadraticFactors2;

	// Call two different functions and compare
	// 1.a X known, search for Y solutions
	C3Vector edge( fieldOfView.GetLowerEdge().GetX(), dummy, sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors );
	// cout << endl;
	// cout << quadraticFactors.aFactor << " " << quadraticFactors.bFactor << " " << quadraticFactors.cFactor 
	//	 << " AND: " << quadraticFactors2.aFactor << " " << quadraticFactors2.bFactor << " " << quadraticFactors2.cFactor << endl;
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);

	// 1.b Y known, search for X solutions
	edge.Set( dummy, fieldOfView.GetLowerEdge().GetY(), sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors );
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);

	// 2.a something else
	coneAxis.Set(coneAxis.GetX() + 1.0, coneAxis.GetY() + 0.5, coneAxis.GetZ());
	edge.Set( dummy, fieldOfView.GetLowerEdge().GetY(), sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors );
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);

	// 2.b something else
	edge.Set( fieldOfView.GetLowerEdge().GetX(), dummy, sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors );
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);

	// 3.a something else
	coneAxis.Set(coneAxis.GetX() + 1.0, coneAxis.GetY() - 2.5, coneAxis.GetZ());
	edge.Set( dummy, fieldOfView.GetLowerEdge().GetY(), sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors );
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_X, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);

	// 3.b something else
	edge.Set( fieldOfView.GetLowerEdge().GetX(), dummy, sliceZ );
	FromConeEquationToQuadraticFactor( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors );
	//	FromConeEquationToQuadraticFactorOld( coneAxis, origin, comptonAngle, Axis_Y, edge, quadraticFactors2 );
	//	assert (quadraticFactors == quadraticFactors2);
}

void
CreateSimpleFOV( CVIPFieldOfView& io_fieldOfView )
{
	int bins_x = 10;
	double xmin = -5.0;
	double xmax =  5.0;
	int bins_y = 10;
	double ymin = -5.0;
	double ymax =  5.0;
	int bins_z = 1;
	double zmin =  8.5;
	double zmax = 11.5;

	io_fieldOfView.SetFieldOfView( bins_x, xmin, xmax, bins_y, ymin, ymax, bins_z, zmin, zmax );
}

// ------------------------------------------------------------------------------------------------------------
// Cone centre is at x = -5, y = 5
// (cone axis is 0, so, ellips is circle, so cone centre = scatter position = origin)
void
CreateSimpleCircle( const double& in_sliceZ, const double& in_circleradius
	, const double& in_origin_x, const double& in_origin_y, C3Vector& io_origin
	, const double& in_axis_x, const double& in_axis_y, C3Vector& io_coneAxis, double& io_comptonAngle)
{
	double tanTh = in_circleradius / in_sliceZ;
	io_comptonAngle = atan( tanTh );

	io_origin.Set(in_origin_x, in_origin_y, 0.0);

	io_coneAxis.Set( in_axis_x, in_axis_y, in_sliceZ );
	double len = io_coneAxis.GetLength();
	double factor = 1.0 / len;
	io_coneAxis = io_coneAxis * factor;

	cout << "CIRCLE, " << io_comptonAngle * 180 / kPI << " degrees: " << " cos(th): " << cos(io_comptonAngle) << endl;
	cout << "origin: " << io_origin << "  ";
	cout << "cone axis: " << io_coneAxis << "  " << endl;
}


// -----------------------------------------------------------------------------------------

void
TestRotateFieldOfView()
{
	cout << endl;
	cout << "*** TestRotateFieldOfView() ***" << endl;

	// TEST CASE..........
	int bins_x = 3;
	double xmin = -1.5;
	double xmax =  1.5;
	int bins_y = 3;
	double ymin = -1.5;
	double ymax =  1.5;
	int bins_z = 1;
	double zmin = -0.5;
	double zmax =  0.5;

	CVIPFieldOfView absoluteFOV;
	absoluteFOV.SetFieldOfView( bins_x, xmin, xmax, bins_y, ymin, ymax, bins_z, zmin, zmax );	

	C3Vector source( -1.0, 1.0, 0.0);
	int absoluteIdx = absoluteFOV.GetVoxelIndex( source );
	cout << "absoluteIdx: " << absoluteIdx << endl;	

	// We rotate the cone, for convenience sake, so we have to rotate the FOV also...
	C3Matrix rotMatrix( Axis_Y, 0.5 * kPI); 	// 90 degrees
	C3Vector rotsource = rotMatrix * source;

	cout << "initial FOV, #bins: " << absoluteFOV.GetNumberOfVoxels() << endl;
	cout << "initial FOV, edges: " << absoluteFOV.GetLowerEdge() << "  " << absoluteFOV.GetUpperEdge() << endl;

	CVIPFieldOfView rotatedFOV( absoluteFOV );
	rotatedFOV.Rotate( rotMatrix );
	unsigned int rotatedIdx = rotatedFOV.GetVoxelIndex( rotsource );

	cout << "rotatedIdx: " << rotatedIdx << endl;		
	cout << "rotated FOV, #bins: " << rotatedFOV.GetNumberOfVoxels() << endl;
	cout << "rotated FOV, edges: " << rotatedFOV.GetLowerEdge() << "  " << rotatedFOV.GetUpperEdge() << endl;

	// With the rotated cone (i.e. the rotated detector hit positions) we solved the cone intersections with the rotated FOV.
	// This gives a certain FOV idx which contains the source voxel.
	// If everything is done right, this FOV idx == rotatedIdx...
	
	// Now we loop over all FOV idx that we found (in this example, there is only one)
	std::vector<unsigned int> fovlist;
	fovlist.push_back( rotatedIdx );
	C3Matrix inverseMatrix( Axis_Y, -0.5 * kPI);
	C3Vector voxelCentre;

	for (int i = 0; i < 1; i++)
	{
		unsigned int fovidx = fovlist[i]; 	// == rotatedIdx
		rotatedFOV.GetVoxelCentre( fovidx, voxelCentre );
	
		voxelCentre = inverseMatrix * voxelCentre;
	
		unsigned int realfovidx = absoluteFOV.GetVoxelIndex( voxelCentre );
//		newfovlist.push_back( realfovidx );
	
		cout << "real idx: " << absoluteIdx << " relative, rotated idx: " 
			 << rotatedIdx << " final, real idx: " << realfovidx << endl;

		assert( realfovidx == absoluteIdx );
	}
}

// -----------------------------------------------------------------------------------------

void TestImageKernel()
{
	cout << "		***  1 ***" << endl;
	cout << endl;
	{
		int size = 3; 	// kernel-size = 3
		int dimension = 2; 	// 2 dimensions
		double sigma = 0.75; 	// In a 3x3 fov array, with voxelsize 1x1 mm^2; this would mean that 95.4% is covered

		CVIPGaussKernel gaussKernel( size, dimension, sigma );
		gaussKernel.PrintGauss();
	}
	cout << endl;
	cout << "		***  2 ***" << endl;
	cout << endl;
	{
		int size = 5; 	// kernel-size = 5
		int dimension = 2; 	// 2 dimensions
		double sigma = 0.75; 	// In a 3x3 fov array, with voxelsize 1x1 mm^2; this would mean that 95.4% is covered

		CVIPGaussKernel gaussKernel( size, dimension, sigma );
		gaussKernel.PrintGauss();
	}

	cout << endl;
	cout << "		***  3 ***" << endl;
	cout << endl;
	// Make a test-image (2D, of size 6x6 pixels, with voxelSize 10)
	{
		CVIPFieldOfView fov;
		int nbins = 6;
		double xmin = -30;
		double xmax = 30;
		fov.SetFieldOfView( nbins, xmin, xmax, nbins, xmin, xmax, 1, 0, 1);
		double voxelSize = (xmax - xmin)/nbins;
		CVIPVector image( fov.GetNumberOfVoxels() );
		double value = 0.0;
		double z = 0.5;
		std::ofstream ofile1("test_oldimage.dat");
		//
		for (int iy = 0; iy < nbins; iy++)
		{
			double y = xmin + (iy + 0.5) * voxelSize;
			for (int ix = 0; ix < nbins; ix++)
			{
				double x = xmin + (ix + 0.5) * voxelSize;
				value = 0.0;
				if ((iy == 1 && ix == 1) || (iy > 2 && ix <= 2))
					value = 1.0;
				else if ((iy == 1 && ix == 4) || (iy == 4 && ix == 4))
					value = 5.0;
				else if ((iy == 1 && ix == 1) || (iy > 2 && ix > 2))
					value = 2.0;
	
				int	index = fov.GetVoxelIndex( ix, iy, 0 );
				image[index] = value;

				// cout << image[index] << "   ";
				ofile1 << x << "    " << y << "     " << z << "     " << image[index] << endl;
			}
			// cout << endl;
		}
		CVIPVector new_image( image );
		// cout << endl;

		int size = 3; 	// kernel-size = 5
		int dimension = 2; 	// 2 dimensions
		double sigma = 0.5; 
		CVIPGaussKernel gaussKernel( size, dimension, sigma );
		gaussKernel.PrintGauss();
		std::ofstream ofile2("test_newimage.dat");
		//
		z = 0.5;
		for (int iy = 0; iy < nbins; iy++)
		{
			double y = xmin + (iy + 0.5) * voxelSize;
			for (int ix = 0; ix < nbins; ix++)
			{
				double x = xmin + (ix + 0.5) * voxelSize;
				int	index = fov.GetVoxelIndex( ix, iy, 0 );
				new_image[index] = gaussKernel.Convolute( image, fov, index );
				// cout << new_image[index] << "   ";
				ofile2 << x << "    " << y << "     " << z << "     " << new_image[index] << endl;
			}
			// cout << endl;
		}
		// cout << endl;		
	}

	// Make a test-image (2D, of size 6x6 pixels, with voxelSize 1)
	{
		CVIPFieldOfView fov;
		int nbinsX = 10;
		int nbinsY = 12;
		double xmin = 0;
		double xmax = 100;
		double ymin = 0;
		double ymax = 120;
		fov.SetFieldOfView( nbinsX, xmin, xmax, nbinsY, ymin, ymax, 1, 0, 1);
		double voxelSizeX = (xmax - xmin)/nbinsX;
		double voxelSizeY = (ymax - ymin)/nbinsY;
		CVIPVector image( fov.GetNumberOfVoxels() );
		double value = 0.0;
		double z = 0.5;
		std::ofstream ofile1("test_oldimage_RS.dat");
		//
		for (int iy = 0; iy < nbinsY; iy++)
		{
			double y = xmin + (iy + 0.5) * voxelSizeY;
			for (int ix = 0; ix < nbinsX; ix++)
			{
				double x = xmin + (ix + 0.5) * voxelSizeX;
				value = 5.0;
				if (iy == 10) {
					if (ix > 2 && ix < 5) value = 10.0;
					if (ix > 5 && ix < 8) value = 10.0;
				}
				else if (iy == 9) { 
					if (ix > 1 && ix < 8) value = 10.0;
				}
				else if (iy == 8) { 
					if (ix == 1 || ix == 7 || ix == 8) value = 10.0;
					if (ix == 2 || ix == 3) value = 0.0;
				}
				else if (iy == 7) { 
					if (ix == 1 || ix == 7 || ix == 8) value = 10.0;
					if (ix == 2) value = 0.0;
				}
				else if (iy >= 4 && iy <= 6) { 
					if (ix == 1 || (ix > 2 && ix < 6) || ix == 7 ) value = 10.0;
				}
				else if (iy == 2 || iy == 3) { 
					if (ix == 0 || (ix > 1 && ix < 5) || ix == 6 ) value = 10.0;
				}
				else if (iy == 1) { 
					if (ix == 1 || ix == 2 || ix == 3 ) value = 10.0;
				}
	
				int	index = fov.GetVoxelIndex( ix, iy, 0 );
				image[index] = value;

				// cout << image[index] << "   ";
				ofile1 << x << "    " << y << "     " << z << "     " << image[index] << endl;
			}
			// cout << endl;
		}
		CVIPVector new_image( image );
		// cout << endl;

		int size = 3; 	// kernel-size = 5
		int dimension = 2; 	// 2 dimensions
		double sigma = 0.75; 
		CVIPGaussKernel gaussKernel( size, dimension, sigma );
		std::ofstream ofile2("test_newimage_RS.dat");
		//
		z = 0.5;
		for (int iy = 0; iy < nbinsY; iy++)
		{
			double y = xmin + (iy + 0.5) * voxelSizeY;
			for (int ix = 0; ix < nbinsX; ix++)
			{
				double x = xmin + (ix + 0.5) * voxelSizeX;
				int	index = fov.GetVoxelIndex( ix, iy, 0 );
				new_image[index] = gaussKernel.Convolute( image, fov, index );
				// cout << new_image[index] << "   ";
				ofile2 << x << "    " << y << "     " << z << "     " << new_image[index] << endl;
			}
			// cout << endl;
		}
		// cout << endl;		
	}

}

// ======================================================================================

void TestVectorsAndMatrices()
{
	// Test unit matrix;
	{
		for (int ival = 1; ival < 4; ival++)
		{
			double diagVal = (double) ival;
			C3Matrix matrix;
			matrix.SetUnity(diagVal);
			for (int icol = 0; icol<3; icol++)
			{
				for (int irow = 0; irow<3; irow++)
				{
					double value = matrix.GetElement(icol, irow);
					if (icol == irow)
					{
						if ( !doubleEquals( value, diagVal, 0.000001 ))
						{
							cout << "unit matrix, icol: " << icol << " irow: " << irow 
								 << " value: " << value << " should be: " << diagVal << endl;
						}
						assert( doubleEquals( value, diagVal, 0.000001 ));
					}
					else
					{
						if ( !doubleEquals( value, 0.0, 0.000001 ))
						{
							cout << "unit matrix, icol: " << icol << " irow: " << irow << " value: " << value << endl;
						}
						assert( doubleEquals( value, 0.0, 0.000001 ));
					}
				}
			}
		}
		cout << "Unit Matrix Test succesful" << endl;
	}
	// Test vector outer product
	{
		C3Vector vec;
		for (int iloop = 0; iloop < 2; iloop++)
		{
			double tgtTrace = 0.0;
			if (iloop == 0)
			{
				vec.Set(1.0, 1.0, 1.0);
				tgtTrace = 3.0;
			}
			else 
			{
				vec.Set(1.0, 2.0, 3.0);
				tgtTrace = 14.0;
			}
			C3Matrix matrix = vec.OuterProduct(vec);
			if ( !matrix.IsSymmetric() )
			{
				cout << "Matrix not symmetric: " << endl;
				cout << matrix << endl;
			}
			assert( matrix.IsSymmetric() );
			double trace = matrix.Trace();
			if ( !doubleEquals( trace, tgtTrace, 0.000001 ) )
			{
				cout << "wrong trace: " << trace << " should be: " << tgtTrace << endl;
			}
			assert (doubleEquals( trace, tgtTrace, 0.000001 ));
	
			for (int icol = 0; icol<3; icol++)
			{
				for (int irow = 0; irow<3; irow++)
				{
					double value = matrix.GetElement(icol, irow);
					if (value != vec.Get(irow) * vec.Get(icol))
					{
						cout << "Wrong outerproduct, icol: " << icol << " vec: " << vec.Get(icol) 
												<< " irow: " << irow << " vec: " << vec.Get(irow) 
												<< " matrix element: " << value << endl;
					}
				}
			}
		}
		cout << "Vector Outer Product Test succesful" << endl;	
	}
}

void TestVipGeometry()
{
	CVIPFieldOfView fieldOfView;
	int bins_x = 11;	double xmin = -5.5;  	double xmax =  5.5;
	int bins_y = 11;	double ymin = -5.5;		double ymax =  5.5;
	int bins_z = 1;		double zmin = -0.5;		double zmax =  0.5;

	fieldOfView.SetFieldOfView( bins_x, xmin, xmax, bins_y, ymin, ymax, bins_z, zmin, zmax );

	// Testing CVipBox
	{
		CVipBox aBox( C3Vector(0.0, 0.0, 0.0), C3Vector(2, 2, 2) );
		assert (  aBox.IsInsideRegion( C3Vector(0.0, 0.0, 0.0) ) );
		assert (  aBox.IsInsideRegion( C3Vector(0.9, 0.9, 0.9) ) );
		assert ( !aBox.IsInsideRegion( C3Vector(1.1, 1.1, 1.1) ) );

		cout << "aBox @:" << aBox.GetPosition() 
			 << " min: " << aBox.GetMinimumPosition() 
			 << " max: " << aBox.GetMaximumPosition() 
			 << endl;
		for (int ivox = 0; ivox < fieldOfView.GetNumberOfVoxels(); ivox++)
		{
			bool isIn =	aBox.IsVoxelInsideGeometry( ivox, fieldOfView );
			if (isIn)
			{
				C3Vector coords;
				fieldOfView.GetVoxelCentre(ivox, coords);
				cout << "        " << ivox << " " << isIn << " @: " << coords << endl;
			}
		}
	}

	// Testing CVipSphere
	{
		CVipSphere aSphere( C3Vector(0.0, 0.0, 0.0), 1.0 ); 
		assert (  aSphere.IsInsideRegion( C3Vector(0.0, 0.0, 0.0) ) );
		assert (  aSphere.IsInsideRegion( C3Vector(0.99, 0.0, 0.0) ) );
		assert (  aSphere.IsInsideRegion( C3Vector(0.0, 0.99, 0.0) ) );
		assert (  aSphere.IsInsideRegion( C3Vector(0.0, 0.0, 0.99) ) );
		assert ( !aSphere.IsInsideRegion( C3Vector(1.01, 0.0, 0.0) ) );
		assert ( !aSphere.IsInsideRegion( C3Vector(0.0, 1.01, 0.0) ) );
		assert ( !aSphere.IsInsideRegion( C3Vector(0.0,  0.0, 1.01) ) );

		cout << "aSphere @:" << aSphere.GetPosition() 
			 << " min: " << aSphere.GetMinimumPosition() 
			 << " max: " << aSphere.GetMaximumPosition() 
			 << endl;
		for (int ivox = 0; ivox < fieldOfView.GetNumberOfVoxels(); ivox++)
		{
			bool isIn =	aSphere.IsVoxelInsideGeometry( ivox, fieldOfView );
			if (isIn)
			{
				C3Vector coords;
				fieldOfView.GetVoxelCentre(ivox, coords);
				cout << "        " << ivox << " " << isIn << " @: " << coords << endl;
			}
		}
	}

	// Testing CVipCircle
	{
		CVipCircle aCircle( C3Vector(0.0, 0.0, 0.0), 1.0 ); 
		assert (  aCircle.IsInsideRegion( C3Vector(0.0, 0.0, 0.0) ) );
		assert (  aCircle.IsInsideRegion( C3Vector(0.99, 0.0, 0.0) ) );
		assert (  aCircle.IsInsideRegion( C3Vector(0.0, 0.99, 0.0) ) );
		assert ( !aCircle.IsInsideRegion( C3Vector(1.01, 0.0, 0.0) ) );
		assert ( !aCircle.IsInsideRegion( C3Vector(0.0, 1.01, 0.0) ) );
		assert ( !aCircle.IsInsideRegion( C3Vector(0.0,  0.0, 0.01) ) );

		cout << "aCircle @:" << aCircle.GetPosition() 
			 << " min: " << aCircle.GetMinimumPosition() 
			 << " max: " << aCircle.GetMaximumPosition() 
			 << endl;
		for (int ivox = 0; ivox < fieldOfView.GetNumberOfVoxels(); ivox++)
		{
			bool isIn =	aCircle.IsVoxelInsideGeometry( ivox, fieldOfView );
			if (isIn)
			{
				C3Vector coords;
				fieldOfView.GetVoxelCentre(ivox, coords);
				cout << "        " << ivox << " " << isIn << " @: " << coords << endl;
			}
		}
	}

	// Testing CVipLine
	{
		CVipLine aLine( C3Vector(0.0, 0.0, 0.0), 4.0 ); 
		assert (  aLine.IsInsideRegion( C3Vector(0.0, 0.0, 0.0) ) );
		assert (  aLine.IsInsideRegion( C3Vector(0.99, 0.0, 0.0) ) );
		assert (  aLine.IsInsideRegion( C3Vector(-0.99, 0.0, 0.0) ) );
		assert ( !aLine.IsInsideRegion( C3Vector(2.01, 0.0, 0.0) ) );

		assert ( !aLine.IsInsideRegion( C3Vector(0.0,  0.01, 0.0) ) );
		assert ( !aLine.IsInsideRegion( C3Vector(0.0, -0.01, 0.0) ) );

		assert ( !aLine.IsInsideRegion( C3Vector(0.0, 0.0,   0.01) ) );
		assert ( !aLine.IsInsideRegion( C3Vector(0.0, 0.0,  -0.01) ) );

		cout << "aLine @:" << aLine.GetPosition() 
			 << " min: " << aLine.GetMinimumPosition() 
			 << " max: " << aLine.GetMaximumPosition() 
			 << endl;
		for (int ivox = 0; ivox < fieldOfView.GetNumberOfVoxels(); ivox++)
		{
			bool isIn =	aLine.IsVoxelInsideGeometry( ivox, fieldOfView );
			if (isIn)
			{
				C3Vector coords;
				fieldOfView.GetVoxelCentre(ivox, coords);
				cout << "        " << ivox << " " << isIn << " @: " << coords << endl;
			}
		}
	}
}


void TestComptonEdge()
{
	double E1_keV, Etot_keV(511.0), angleRad, angleDeg;
	for (int iA = 0; iA < 360; iA++)
	{	
		angleDeg = (double) iA;
		angleRad = angleDeg * kPI/180.0;
		
		bool okA = VIPUtils::GetComptonEnergyFromAngleRad(angleRad, Etot_keV, E1_keV );	
		assert( okA );
		if (okA)
			cout << "Valid angle: " << angleDeg << " gives E1: " << E1_keV << endl;
		else
			cout << "Invalid angle: " << angleDeg << " gives E1: " << E1_keV << endl;
	
		bool okE = VIPUtils::GetComptonAngleDegrees( E1_keV, Etot_keV, angleDeg );
		if (okE)
			cout << "Valid E1: " << E1_keV << " gives angle: " << angleDeg << endl;
		else
			cout << "Invalid E1: " << E1_keV << " gives angle: " << angleDeg << endl;
	}
}

void TestSolveSetOfEquations_grade3()
{
    int ndim(4);
	CVIPMatrix origA(ndim, ndim);
	CVIPVector origb(ndim);

	// TEST 1, DONE ON PAPER:
	/*
	E =  100, 200, 500, 1000
	E' = 111.01, 214.08, 536.25, 1120
	p = 10, 1, 10e-4, 10e-8
	*/
	origb[0] = 111.01;
	origb[1] = 214.08;
	origb[2] = 536.25;
	origb[3] = 1120.0;

	// Do Gaussian Elimination
	GaussianElimination gaussianElimination(ndim);

	TH1D* h1 = new TH1D("h1", "", 1400, 0, 1400);
	h1->SetMinimum(0.0);
	h1->SetMaximum(1400.0);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TF1* fitPol2 = ((TF1 *)(gROOT->GetFunction("pol2")));

    double xarr[ndim]; 
	xarr[0] = 100;
	xarr[1] = 200; 
	xarr[2] = 500;
	xarr[3] = 1000;
    double coeff;
	double p0, p1;

	for (int idim = 0; idim < ndim; idim++)
	{
		double x = xarr[idim];
		for (int icolumn = 0; icolumn < ndim; icolumn++)
		{
			double coeff = 1;
			for (int ic = 0; ic < icolumn; ic++)
				coeff = coeff*x;
			origA.SetElement(idim, icolumn, coeff );
		}
	}
	// Do Gaussian Elimination
	gaussianElimination.SetVectorB( origb );
	gaussianElimination.SetMatrixA( origA );
	gaussianElimination.PrintEquations();
	//	int fck; cin >> fck;
	gaussianElimination.SolveSet( false );

	CVIPVector solutionP(ndim);
	gaussianElimination.GetSolutions( solutionP );

	cout << "SOLUTION P: " << solutionP << endl;
	for (int irow = 0; irow < ndim; irow++)
	{
		double y_axis(0.0);
		for (int icol = 0; icol < ndim; icol++)
		{
			y_axis += solutionP[icol] * origA.GetElement(irow, icol);
		}
		cout << irow << " gives y= " << y_axis << " (real-value: " << origb[irow] << ")" << endl;
	}

	// =======================================================

	origb[0] = m_peakFitBa133_303;
	origb[1] = m_peakFitBa133_356;
	origb[2] = m_peakFitNa22_511;
	origb[3] = m_peakFitNa22_1274;

	int iset(0);
	for (int ievent = 0; ievent < 2; ievent++)
	{
		if (ievent == 0)
		{
			xarr[0] = 303.807; xarr[1] = 356.499; xarr[2] = 514.347; xarr[3] = 1206.79;
		}
		else if (ievent == 1)
		{
			xarr[0] = 304.301; xarr[1] = 357.273; xarr[2] = 513.884; xarr[3] = 1186.26;
		}

		for (int idim = 0; idim < ndim; idim++)
		{
			double x = xarr[idim];
			for (int icolumn = 0; icolumn < ndim; icolumn++)
			{
				double coeff = 1;
				for (int ic = 0; ic < icolumn; ic++)
					coeff = coeff*x;
				origA.SetElement(idim, icolumn, coeff );
			}
		}
		// Do Gaussian Elimination
		gaussianElimination.SetVectorB( origb );
		gaussianElimination.SetMatrixA( origA );
		gaussianElimination.PrintEquations();
		//	int fck; cin >> fck;
		gaussianElimination.SolveSet( false );

		gaussianElimination.GetSolutions( solutionP );

		cout << "SOLUTION P: " << solutionP << endl;
		double gr_x[4], gr_y[4], gr_delX[4], gr_delY[4];
		for (int irow = 0; irow < ndim; irow++)
		{
			double y_axis(0.0);
			for (int icol = 0; icol < ndim; icol++)
			{
				y_axis += solutionP[icol] * origA.GetElement(irow, icol);
			}
			cout << irow << " gives y= " << y_axis << " (real-value: " << origb[irow] << ")" << endl;
			gr_x[irow] = origA.GetElement(irow, 1);
			gr_delX[irow] = 0.025*gr_x[irow];
			gr_y[irow] = y_axis;
			gr_delY[irow] = 0.025*gr_y[irow];
		}

		TCanvas* m_1 = new TCanvas("m_1", "", 700, 700);
		h1->Draw();

		TGraphErrors* graph1 = new TGraphErrors(ndim, gr_x, gr_y, gr_delX, gr_delY);
		// graph1->Draw("ALP");
		graph1->Draw("SAME LP");

		TF1* myFunc = new TF1("fitFunc", VIPUtils::PolynomialGrade3, 0, 1400, 4);
		myFunc->SetParameters(solutionP[0], solutionP[1], solutionP[2], solutionP[3]);
		myFunc->SetLineColor(kMagenta);
		myFunc->Draw("SAME");

		CalculateP0_P1_Approx(gr_x, gr_y, p0, p1);
		fitPol2->SetParameter( 0, p0 );
		fitPol2->SetParameter( 1, p1 );
		fitPol2->SetParameter( 2, 0.00001 );

		fitPol2->SetParLimits( 0, -100, 100);
		fitPol2->SetParLimits( 1, 0.0, 1.5);
		fitPol2->SetParLimits( 2, 0.0, 0.001);

		graph1->Fit("pol2", "Q", "", 0, 1300);
		//	TF1* fitPol2 = graphE_2->GetFunction("pol2");
		//	cout << "P0: " << fitPol2->GetParameter(0) << endl;
		fitPol2->SetLineColor(kRed);

		gStyle->SetStatX(0.5);
		gStyle->SetStatY(0.8);

		TString canvasname("canvas_");
		canvasname += iset;
		canvasname += ".png";
		m_1->Print(canvasname, "png");

		delete m_1;
		delete graph1;

		iset++;
	}
}

// ===================================================================

void
CalculateP0_P1_Approx(const double* in_x, const double* in_y, double& out_p0, double& out_p1)
{
	double delY = in_y[1] - in_y[0];
	double delX = in_x[1] - in_x[0];
	out_p1 = delY/delX;
	out_p0 = in_y[0] - in_x[0]*out_p1;
	//	cout << "delX: " << delX << " delY: " << delY << " p0: " << out_p0 << " p1: " << out_p1 << endl;
}

// ===================================================================

void TestSolveSetOfEquations_grade4()
{
    int ndim(5);
	CVIPMatrix origA(ndim, ndim);
	CVIPVector origb(ndim);

	// TEST 1, DONE ON PAPER:
	/*
	E =  10,    100,      200,      500,      1000
	E' = 20.01, 111.0101, 214.0816, 536.3125, 1121
	p =  10, 1, 10e-4, 10e-8, 10e-12
	*/
	origb[0] =  20.01 + 1e-5 + 1e-8;
	origb[1] = 111.0101;
	origb[2] = 214.0816;
	origb[3] = 536.3125;
	origb[4] = 1121.0;

	// Do Gaussian Elimination
	GaussianElimination gaussianElimination(ndim);

	TH1D* h1 = new TH1D("h1", "", 1400, 0, 1400);
	h1->SetMinimum(0.0);
	h1->SetMaximum(1400.0);

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

    double xarr[ndim];
	xarr[0] = 10;
	xarr[1] = 100;
	xarr[2] = 200;
	xarr[3] = 500;
	xarr[4] = 1000;
    double coeff;
	double p0, p1;

	for (int idim = 0; idim < ndim; idim++)
	{
		double x = xarr[idim];
		for (int icolumn = 0; icolumn < ndim; icolumn++)
		{
			double coeff = 1;
			for (int ic = 0; ic < icolumn; ic++)
				coeff = coeff*x;
			origA.SetElement(idim, icolumn, coeff );
		}
	}
	// Do Gaussian Elimination
	gaussianElimination.SetVectorB( origb );
	gaussianElimination.SetMatrixA( origA );
	gaussianElimination.PrintEquations();
	//	int fck; cin >> fck;
	gaussianElimination.SolveSet( false );

	CVIPVector solutionP(ndim);
	gaussianElimination.GetSolutions( solutionP );

	cout << "SOLUTION P: " << solutionP << endl;
	for (int irow = 0; irow < ndim; irow++)
	{
		double y_axis(0.0);
		for (int icol = 0; icol < ndim; icol++)
		{
			y_axis += solutionP[icol] * origA.GetElement(irow, icol);
		}
		cout << irow << " gives y= " << y_axis << " (real-value: " << origb[irow] << ")" << endl;
	}
}

// =============================================

void TestVectorSubtraction()
{
	C3Vector position1( 10, 10, 1 );
	C3Vector position2(  5, 1, 0.5 );

	C3Vector deltaPos = position2 - position1;
	double distance1 = deltaPos.GetLength();
	double distance2 = (position2-position1).GetLength();
	cout << "should be equal: " << distance1 << " and " << distance2 << " delta: " 
		  << " also: " << deltaPos << " and " << (position2-position1) << endl;

	deltaPos = position1 - position2;
	distance1 = deltaPos.GetLength();
	distance2 = (position1-position2).GetLength();
	cout << "should be equal: " << distance1 << " and " << distance2 << " delta: " 
		  << " also: " << deltaPos << " and " << (position1-position2) << endl;
}

// =============================================

void TestVIPRandom()
{
    TH1D* h1 = new TH1D("h1", "uniform", 100, 0., 2.0);
    TH1D* h2 = new TH1D("h2", "gauss", 100, 0., 100.0);

	double umin(0.8), umax(1.7);
	cout << "UNIFORM: " << umin << " " << umax << endl;

    CVIPRandomUniform uniformDist(umin, umax);
    for (int i = 0; i < 1000; i++)
    {
        double value = uniformDist.GetNewValue();
        h1->Fill(value);
    }

	double gmean(54.0), gsigma(8.0);
	cout << "GAUSS: " << gmean << " " << gsigma << endl;

    int seed(1);
    CVIPRandomGauss gaussDist(gmean, gsigma, seed);
    for (int i = 0; i < 10000; i++)
    {
        double value = gaussDist.GetNewValue();
        h2->Fill(value);
    }

    TCanvas* m_1 = new TCanvas("m_1", " ", 1000, 500);
    m_1->Divide(2, 1);
    m_1->cd(1);
    h1->Draw();
    m_1->cd(2);
    gStyle->SetOptFit(1);   // moron!
    h2->Draw();
    h2->Fit("gaus");
    m_1->Print("myrandoms.png", "png");

}

void AuxGetTheRandoms( double& aUniform, double& aGauss)
{
	CVIPRandomUniform uniformDist(2011.7, 2080.1);
	aUniform = uniformDist.GetNewValue();

	CVIPRandomGauss gaussDist(542.0, 16.0);
	aGauss = gaussDist.GetNewValue();
}

void TestVIPRandom2()
{
    TH1D* h1_b = new TH1D("h1_b", "uniform", 100, 2000.0, 2100.0);
    TH1D* h2_b = new TH1D("h2v", "gauss", 100, 450.0, 650.0);

	double umin(2011.7), umax(2080.1);
	cout << "UNIFORM: " << umin << " " << umax << endl;

	double gmean(542.0), gsigma(16.0);
	cout << "GAUSS: " << gmean << " " << gsigma << endl;

	double aFlat, aGauss;
    for (int i = 0; i < 10000; i++)
    {
		AuxGetTheRandoms( aFlat, aGauss);
        h1_b->Fill(aFlat);
        h2_b->Fill(aGauss);
	}
    TCanvas* m_2 = new TCanvas("m_2", " ", 1000, 500);
    m_2->Divide(2, 1);
    m_2->cd(1);
    h1_b->Draw();
    m_2->cd(2);
    gStyle->SetOptFit(1);   // moron!
    h2_b->Draw();
    h2_b->Fit("gaus");
    m_2->Print("myrandoms_2.png", "png");

}



