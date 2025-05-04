#include "CVIPFieldOfView.h"
#include "CVIPImage.h"
  
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

const double mu_value_material = 0.00968; // mm-1

double my_mu_value(0.0);
string unitstr("mm^-1");

/*
 * EXPLANATION:
 *
 * IWDM article - 2016:
 * ====================
 *
 * we used a basic attenuation correction procedure, where all pixels inside the phantom are assigned an attenuation correction factor of 0.0096 mm^-1, and the weights of the LORs are normalized in proportion to how much of their length is traversing the phantom.
 *
 * -----------------------------------------------------
 *
 * Scatter and attenuation corrections for a PEM detector using List-Mode OSEM  (Cláudia S. Ferreira, et al.)
 * ============================================================================
 * https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6829264
 *
 * Beer-Lambert law: I/I_0 = exp(-1.0 * μ * x) ==> ACF = exp(μ * x)
 *
 * where μ is the linear attenuation coefficient, which for 511 keV photons in breast tissue is ~ 0.0958 cm^-1 .
 *  (Machiel: = 0.00958 mm^-1)
 *
 * ACF = Attenuation correction factor
 *
 * Formula 3 (LM-OSEM):
 * x(n) = x(n) SIGMA_event [ ACF * tij/ sigma_all_FOV_bins_for_this_event x(FOV_bin)]
 *     (event = LOR)
 *
 * NOTE: ACF = exp(1.0 * mu * x)  (NOTE the plus sign!!!)
 *
 * ------------------------------------------------------
 *
 * CC_OSEM:
 * =========
 * double att = m_attMap.GetVoxelValue( index );
 * double attcor = exp(-1.0 * att * distance);
 *
 * NOTE: attcor = exp(-1.0 * att * distance)  (NOTE the minus sign!!!)
 *
 * sensitivity_j *= attcor;
 *
 * lambda_j^{l+1} = (lambda^{l} / s_j)
 *     * SUM_i^{Nevts} [ (delta_{ji} * T_{ij}) /  (SUM_{k=1}^{M_FOVbinsOnThisEvt_i} delta_{ki} T_{ik} lambda_k) ]
*/

void FromOSEMImage( const CVIPFieldOfView& fieldOfView );
void FromGeometry( const CVIPFieldOfView& fieldOfView  );
bool IsInSphere( const double& xpos, const double& ypos, const double& zpos, const double& sphere_R );
bool IsInDisk( const double& xpos, const double& ypos, const double& zpos, const double& sphere_R
    , const double& disk_zmin, const double& disk_zmax );
bool IsInBox( const double& xpos, const double& ypos, const double& zpos
             , const double& box_xmin, const double& box_xmax
             , const double& box_ymin, const double& box_ymax
             , const double& box_zmin, const double& box_zmax );

// =======================================================================

int main( int argc, char *argv[] )
{
    my_mu_value = mu_value_material;
    int ichoice(0);
    while (ichoice != 1 && ichoice != 2)
    {
        cout << "make attenuation map from geometry (1) or from image file (2)" << endl;
        cin >> ichoice;
    }

	CVIPFieldOfView fieldOfView;
    try
    {
        fieldOfView.Initialize( "osem_fov_parameters.conf" );
    }
    catch (int e)
    {
        cout << "An exception occurred. Exception Nr. " << e << endl;
        return 1;
    }
    
    if (ichoice == 1)
        FromGeometry(fieldOfView);
    else if (ichoice == 2)
        FromOSEMImage(fieldOfView); 
}

// =======================================================

void FromGeometry( const CVIPFieldOfView& in_fieldOfView )
{
    enum GEOM_TYPE
    {
        GEOM_NONE = 0,
        GEOM_SPHERE = 1,
        GEOM_BOX = 2,
        GEOM_DISK = 3, 
        GEOM_SPHERE_IN_SPHERE = 4,
        GEOM_HALF_SPHERE = 5, 
        GEOM_HALF_SPHERE_IN_HALF_SPHERE = 6,
        GEOM_DISK_IN_DISK = 7 
    };
    const int MAX_GEOM_TYPE = GEOM_DISK_IN_DISK;
    
    enum ORIENTATION
    {
        ORIENT_NONE = 0,
        ORIENT_PLUS_X,
        ORIENT_MIN_X,
        ORIENT_PLUS_Y,
        ORIENT_MIN_Y,
        ORIENT_PLUS_Z,
        ORIENT_MIN_Z
    };
    ORIENTATION orientation(ORIENT_NONE);

    double sphere_R, sphere_R_in;
    double box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax;
    double my_mu_value_out(0.0), my_mu_value_in(0.0);

    cout << "Select if the source object is a sphere or a box." << endl;
    cout << " Sphere = 1, Box = 2, Disk = 3, Sphere in Sphere = 4" << endl;
    cout << " Half-Sphere = 5, Half-Sphere in Half-Sphere = 6" << endl;
    cout << " Disk in Disk = 7" << endl;
    int geomType; cin >> geomType;
    if (geomType > MAX_GEOM_TYPE) geomType = GEOM_NONE;

    if (geomType == GEOM_SPHERE)
    {
        cout << "Give the radius of the sphere (mm)" << endl;
        cin >> sphere_R;
    }
    else if (geomType == GEOM_BOX)
    {
        cout << "Give the size of the box (mm)" << endl;
        cout << "Minimum X and maximum X: " << endl; cin >> box_xmin >> box_xmax;
        cout << "Minimum Y and maximum Y: " << endl; cin >> box_ymin >> box_ymax;
        cout << "Minimum Z and maximum Z: " << endl; cin >> box_zmin >> box_zmax;
    }
    else if (geomType == GEOM_DISK)
    {
        cout << "Give the radius of the disk (mm)" << endl;
        cin >> sphere_R;
        cout << "Minimum Z and maximum Z (mm): " << endl; cin >> box_zmin >> box_zmax;
    }
    else if (geomType == GEOM_SPHERE_IN_SPHERE)
    {
        cout << "Give the radius of the outer sphere (mm)" << endl;
        cin >> sphere_R;
        cout << "Give the radius of the inner sphere (mm)" << endl;
        cin >> sphere_R_in;
        cout << "give value for 'mu' for outer sphere (mm^-1)" << endl;
        cin >> my_mu_value_out;
        cout << "give value for 'mu' for inner sphere (mm^-1)" << endl;
        cin >> my_mu_value_in;
    }
    else if (geomType == GEOM_HALF_SPHERE)
    {
        cout << "Give the radius of the half sphere (mm)" << endl;
        cin >> sphere_R;
        cout << "Give orientation of half sphere: 1 = +X, 2 = -X, 3 = +Y, 4 = -Y" << endl;
        int tmp; cin >> tmp;
        if (tmp <= ORIENT_MIN_Y)
            orientation = static_cast<ORIENTATION> (tmp);
    }
    else if (geomType == GEOM_HALF_SPHERE_IN_HALF_SPHERE)
    {
        cout << "Give the radius of the outer half sphere (mm)" << endl;
        cin >> sphere_R;
        cout << "Give the radius of the inner half sphere (mm)" << endl;
        cin >> sphere_R_in;
        cout << "Give orientation of half sphere: " << endl;
        cout << "1 = +X, 2 = -X, 3 = +Y, 4 = -Y, 5 = +Z, 6 = -Z" << endl;
        int tmp; cin >> tmp;
        if (tmp <= ORIENT_MIN_Z)
            orientation = static_cast<ORIENTATION> (tmp);
cout << "CHECK, orientation = " << orientation << endl;
        
        cout << "give value for 'mu' for outer sphere (mm^-1)" << endl;
        cin >> my_mu_value_out;
        cout << "give value for 'mu' for inner sphere (mm^-1)" << endl;
        cin >> my_mu_value_in;
    }
	//
    else if (geomType == GEOM_DISK_IN_DISK)		// GEOM_DISK_IN_DISK = 7 
    {
        cout << "Give the radius of the outer disk (mm)" << endl;
        cin >> sphere_R;
        cout << "Give the radius of the inner disk (mm)" << endl;
        cin >> sphere_R_in;
        cout << "give value for 'mu' for outer disk (mm^-1)" << endl;
        cin >> my_mu_value_out;
        cout << "give value for 'mu' for inner disk (mm^-1)" << endl;
        cin >> my_mu_value_in;
        cout << "Minimum Z and maximum Z (mm): " << endl; cin >> box_zmin >> box_zmax;
	}
    else
    {
        cout << "wrong geomType: " << geomType << endl;
        exit(1);
    }
    
    // Mu value
    //
    if (   geomType != GEOM_SPHERE_IN_SPHERE 
		&& geomType != GEOM_HALF_SPHERE_IN_HALF_SPHERE
		&& geomType != GEOM_DISK_IN_DISK )
    {
        cout << "using mu value (water): " << my_mu_value << " " << unitstr << endl;
        cout << "give different if you want to use a different mu, '0' = using current value" << endl;
        int tmp; cin >> tmp;
        if (tmp > 0)
            my_mu_value = tmp;

		cout << "my_mu_value: " << my_mu_value << endl;
		int isok(0); cout << "ok? CTRL-C if something is wrong. Hit any key else." << endl; cin >> isok;
    }
	else
	{
        cout << "value for 'mu' for outer sphere (mm^-1): " << my_mu_value_out << endl;
        cout << "value for 'mu' for inner sphere (mm^-1): " << my_mu_value_in << endl;
	}
    
    double shiftInZ(0.0);
    cout << "Give the shift in Z of the sphere/half-sphere/whatever 0.0 = none" << endl;
    cin >> shiftInZ;

    int nvoxels = in_fieldOfView.GetNumberOfVoxels();
    cout << "nvoxels: " << nvoxels << endl;
    cout << "nX, nY, nZ: " << in_fieldOfView.GetNumVoxelsX() << " "
                           << in_fieldOfView.GetNumVoxelsY() << " "
                           << in_fieldOfView.GetNumVoxelsZ() << endl;

    CVipImage image( nvoxels );

    int ix, iy, iz;
    int numX = in_fieldOfView.GetNumVoxelsX();
    int numY = in_fieldOfView.GetNumVoxelsY();
    for (int iv = 0; iv < nvoxels; iv++)
    {
        in_fieldOfView.GetVoxelIndices( iv, ix, iy, iz);
        double xpos(0.0), ypos(0.0), zpos(0.0);
        bool isValid = in_fieldOfView.GetVoxelCentre(ix, iy, iz, xpos, ypos, zpos);
        bool insideSource(false);
        if (isValid)
        {
            if ( abs(shiftInZ) > 0.0 )
                zpos = zpos - shiftInZ;
            
            if (geomType == GEOM_SPHERE)
                insideSource = IsInSphere( xpos, ypos, zpos, sphere_R );
            else if (geomType == GEOM_BOX)
                insideSource = IsInBox( xpos, ypos, zpos
                    , box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax );
            else if (geomType == GEOM_DISK)
                insideSource = IsInDisk( xpos, ypos, zpos, sphere_R, box_zmin, box_zmax );
			//
            else if (geomType == GEOM_DISK_IN_DISK)
			{
                insideSource = IsInDisk( xpos, ypos, zpos, sphere_R_in, box_zmin, box_zmax ); 
                if ( insideSource) 
                    my_mu_value = my_mu_value_in;
                else 
				{
					insideSource = IsInDisk( xpos, ypos, zpos, sphere_R, box_zmin, box_zmax );
                	if ( insideSource) 
                    	my_mu_value = my_mu_value_out;
				}
			}
   			// 
            else if (geomType == GEOM_SPHERE_IN_SPHERE)
            {
                insideSource = IsInSphere( xpos, ypos, zpos, sphere_R_in );
                if (insideSource)
                {
                    my_mu_value = my_mu_value_in;
                }
                else 
                {
                    insideSource = IsInSphere( xpos, ypos, zpos, sphere_R );
                    if (insideSource) 
                    {
                        my_mu_value = my_mu_value_out;
                    }
                }
            }
            
            else if (geomType == GEOM_HALF_SPHERE)
            {
                insideSource = IsInSphere( xpos, ypos, zpos, sphere_R );
                insideSource = insideSource && (orientation != ORIENT_PLUS_X || xpos > 0);
                insideSource = insideSource && (orientation != ORIENT_MIN_X || xpos < 0);
                insideSource = insideSource && (orientation != ORIENT_PLUS_Y || ypos > 0);
                insideSource = insideSource && (orientation != ORIENT_MIN_Y || ypos < 0);
                insideSource = insideSource && (orientation != ORIENT_PLUS_Z || zpos > 0);
                insideSource = insideSource && (orientation != ORIENT_MIN_Z || zpos < 0);
            }
            else if (geomType == GEOM_HALF_SPHERE_IN_HALF_SPHERE)
            {
                insideSource = IsInSphere( xpos, ypos, zpos, sphere_R_in );
                if (insideSource)
                {
                    my_mu_value = my_mu_value_in;
                }
                else 
                {
                    insideSource = IsInSphere( xpos, ypos, zpos, sphere_R );
                    if (insideSource) 
                    {
                        my_mu_value = my_mu_value_out;
                    }
                }
                if (insideSource)
                {
                    insideSource = insideSource && (orientation != ORIENT_PLUS_X || xpos > 0);
                    insideSource = insideSource && (orientation != ORIENT_MIN_X || xpos < 0);
                    insideSource = insideSource && (orientation != ORIENT_PLUS_Y || ypos > 0);
                    insideSource = insideSource && (orientation != ORIENT_MIN_Y || ypos < 0);
                    insideSource = insideSource && (orientation != ORIENT_PLUS_Z || zpos > 0);
                    insideSource = insideSource && (orientation != ORIENT_MIN_Z || zpos < 0);
                }
            }
        }

        if (!insideSource)
            image.SetVoxelValue( iv, 0.0 );
        else
        {
            // cout << "filling: " << iv << ": " << ix << ":" << iy << ":" << iz
            //      << " @: " << xpos << " " << ypos << " " << zpos 
            //      << " mu: " << my_mu_value << endl;
            
            image.SetVoxelValue( iv, my_mu_value );
        }
    }

    string fname_new( "ATTENUATION_map.img_NEW" );
    FILE_FORMAT fileformat( FFORMAT_BINARY );
    image.Write( fname_new, in_fieldOfView, fileformat );
}

// =======================================================

bool IsInSphere( const double& xpos, const double& ypos, const double& zpos, const double& sphere_R )
{
    double distance = sqrt( xpos*xpos + ypos*ypos + zpos*zpos );
    bool isinside = distance < sphere_R;
    return isinside;
}

// =======================================================

bool IsInDisk( const double& xpos, const double& ypos, const double& zpos, const double& sphere_R
    , const double& disk_zmin, const double& disk_zmax )
{
    double distance = sqrt( xpos*xpos + ypos*ypos );
    bool isinside = distance < sphere_R;
    if (isinside)
    {
        isinside = disk_zmin < zpos && zpos < disk_zmax;
    }

    return isinside;
}

// =======================================================

bool IsInBox( const double& xpos, const double& ypos, const double& zpos
             , const double& box_xmin, const double& box_xmax
             , const double& box_ymin, const double& box_ymax
             , const double& box_zmin, const double& box_zmax )
{
    return (    box_xmin < xpos && xpos < box_xmax
             && box_ymin < ypos && ypos < box_ymax
             && box_zmin < zpos && zpos < box_zmax );
}

// =======================================================

// Old code where every pixel of an input image is converted to a "water-attenuated pixel"
//  (and hence, given the value "mu_value_material")
//
void FromOSEMImage( const CVIPFieldOfView& in_fieldOfView )
{
	cout << "WARNING WARNING WARNING!!!!" << endl;
	cout << "some bins in the original image might be 0, this will result in holes in the attenuation map, with severe consequences!!!!" << endl;
	cout << "(end of WARNING)" << endl;

	string fname;
	cout << "give starting image (*.img) filename" << endl;
	cin >> fname;
    
    cout << "give value for 'mu' (" << unitstr << ")" << endl;
    cin >> my_mu_value;    

	CVipImage image;
	FILE_FORMAT fileformat( FFORMAT_BINARY );
	DATA_FORMAT dataformat( DFORMAT_FLOAT );
	image.Read( fname, in_fieldOfView, fileformat, dataformat);

	int nvoxels = image.GetSize();

	int nvoxels_check = in_fieldOfView.GetNumberOfVoxels();
	if ( nvoxels != nvoxels_check )
	{
		cout << "Image size: " << nvoxels << endl;
		cout << "FOV size: " << nvoxels_check << endl;
	}
	assert( nvoxels == nvoxels_check );

	cout << "nvoxels: " << nvoxels << endl;
	cout << "give minimum signal threshold" << endl;
	double threshold;
	cin >> threshold;

	int ix, iy, iz;
	int numX = in_fieldOfView.GetNumVoxelsX();
	int numY = in_fieldOfView.GetNumVoxelsY();
	for (int iv = 0; iv < nvoxels; iv++) 
	{
		in_fieldOfView.GetVoxelIndices( iv, ix, iy, iz);

		double val = image.GetVoxelValue( iv );
		if (val < threshold) // || ix < 2 || ix > numX-3 || iy < 2 || iy > numY-3)
			image.SetVoxelValue( iv, 0.0 );
		else
        {
			//   image.SetVoxelValue( iv, mu_value_material );
            image.SetVoxelValue( iv, my_mu_value );
        }
	}

	string fname_new( "ATTENUATION_map.img_NEW" );
	image.Write( fname_new, in_fieldOfView, fileformat );
}



