#pragma once

#ifndef _VIPCOMMON_UTILS_H__
#define _VIPCOMMON_UTILS_H__

#include "VIPconstants.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class C3Vector;

enum AXIS
{
	Axis_X = 0,
	Axis_Y,
	Axis_Z
};

void PrintAxis( AXIS in_axis );

bool doubleEquals(double left, double right, double epsilon);

bool fileExists(const std::string& filename);

// ---------------- class QuadraticFactors -----------------

class QuadraticFactors
{
public:
	// ctor and dtor
					QuadraticFactors() : m_aFactor(0), m_bFactor(0), m_cFactor(0) {}
					QuadraticFactors(const double& a, const double& b, const double& c)
						: m_aFactor(a), m_bFactor(b), m_cFactor(c) {}
					~QuadraticFactors() {}

	// get functions
	inline double	GetaFactor() const { return m_aFactor; }
	inline double	GetbFactor() const { return m_bFactor; }
	inline double	GetcFactor() const { return m_cFactor; }

	inline void		SetaFactor( const double& in_aFactor ) { m_aFactor = in_aFactor; }
	inline void		SetbFactor( const double& in_bFactor ) { m_bFactor = in_bFactor; }
	inline void		SetcFactor( const double& in_cFactor ) { m_cFactor = in_cFactor; }

	inline bool		operator==(const QuadraticFactors& in_obj) const 		// ==
	{
		return (   doubleEquals( m_aFactor, in_obj.GetaFactor(), 0.000001 )
				&& doubleEquals( m_bFactor, in_obj.GetbFactor(), 0.000001 )
				&& doubleEquals( m_cFactor, in_obj.GetcFactor(), 0.000001 ) );
	}

private:

	// Prevent use of copy and assign constructors
						QuadraticFactors(const QuadraticFactors& in_obj);
	QuadraticFactors& 	operator= (const QuadraticFactors& in_obj);

	// data
	double 			m_aFactor;
	double 			m_bFactor;
	double 			m_cFactor;
};

// TODO, put next functions inside a namespace?
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
void WriteString(std::ostream& io_os, const std::string& in_string);
void ReadString(std::istream& in_is, std::string& io_string);

void bubbleSort( int size, double* x );

int
solveQuadraticEquation(const QuadraticFactors& in_quadraticFactors, double& out_x1, double& out_x2);

double
GetPhiFromXandY(const double& in_X, const double& in_Y);

bool
RayBoxIntersection(const C3Vector& in_origin, const C3Vector& in_direction
				 , const C3Vector& lowerEdge, const C3Vector& upperEdge
				 , double& io_factor_enter, double& io_factor_exit );

double
NextRayIntersection(  const C3Vector& in_origin, const C3Vector& in_direction
					, const C3Vector& in_minBounds, const C3Vector& in_maxBounds );

bool
IsInsideBounds(const C3Vector& in_minBounds, const C3Vector& in_maxBounds, const C3Vector& in_position);


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

namespace VIPUtils
{
    // Use & reference to avoid copy constructor for ifstream (and because it's better anyway)
    template<typename T> inline void DoReadBinary(std::ifstream& infile,  T& in_tmp)
    {
        infile.read((char*) &in_tmp, sizeof(in_tmp));
    }
    
    // Use & reference to avoid copy constructor for ifstream (and because it's better anyway)
    template<typename T> inline void DoWriteBinary(std::ofstream& outfile, const T& in_tmp)
    {
        outfile.write((char*) &in_tmp, sizeof(in_tmp));
    }
    
	void	GetRandomSolidAngle(double& io_phi, double& io_theta, bool doFourPi);

	void	GetRandomPositionOnSolidAngle( const double& in_phi, const double& in_theta
					, const double& in_zMin, const double& in_zSize, const C3Vector& in_originPosition
					, C3Vector& out_randomPosition );

	double	GetPseudoSolidAngleProbability( const C3Vector& in_originCoords, const C3Vector& in_surfaceCoords );

	double	Gauss( const C3Vector& in_mean, const double& in_sigma, const C3Vector& in_coords, int in_dim );

    double  minAbsDifferenceTwoAnglesRad( const double& in_angle1, const double& in_angle2 );
	bool 	GetGeomAngleDegrees(const C3Vector& v1, const C3Vector& v2, double& out_angleGeomDegrees);

	bool	GetPositionAngle2D( const double& in_x, const double& in_y, double& out_angleGeom);

	double	GetDistance(const C3Vector& v1, const C3Vector& v2);

	// WARNING!!! THIS FUNCTION WORKS WITH keV!!!!!
	bool	GetComptonAngleDegrees(const double& E1_keV, const double& Etot_keV
                        , double& out_angleComptonDegrees );
	bool	GetComptonEnergyFromAngleRad(const double& in_angleComptonRadians
					, const double& Etot_keV, double& out_E1_keV );

	// parse COMMA separated substrings
	void	GetFileNames ( const std::string& in_string, std::vector<std::string>& out_fileNames ); 

	double PolynomialGrade3( double* in_x, double* in_par );
	double PolynomialGrade4( double* in_x, double* in_par );
    
    //
    // The next function calculate the LOR between two hits X1 and X2 
    //  and then the distance this LOR is traversing through a sphere (like a source)
    //  where the sphere has a radius and a centre position. 
    //  It returns true if it found 2 intersection points of the LOR with the sphere and 
    //  the return variable "out_distance" will be the distance of the LOR that lies within the sphere
    //
    // This can be used to calculate attenuation correction for a spherical source.
    // 
    bool 
    CalculateIntersectionDistance(const C3Vector& in_X1, const C3Vector& in_X2
        , const C3Vector in_Centre, const double& in_radius, double& out_distance);
};


#endif




