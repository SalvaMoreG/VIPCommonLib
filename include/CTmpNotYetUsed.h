
/*
void
Bresenham3DLinePlotting( const C3Vector& in_P1, const C3Vector& in_P2, int in_precision )
{

%  Generate X Y Z coordinates of a 3D Bresenham's line between
%  two given points.
%
%  A very useful application of this algorithm can be found in the
%  implementation of Fischer's Bresenham interpolation method in my
%  another program that can rotate three dimensional image volume
%  with an affine matrix:
%  http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=21080
%
%  Usage: [X Y Z] = bresenham_line3d(P1, P2, [precision]);
%
%  P1	- vector for Point1, where P1 = [x1 y1 z1]
%
%  P2	- vector for Point2, where P2 = [x2 y2 z2]
%
%  precision (optional) - Although according to Bresenham's line
%	algorithm, point coordinates x1 y1 z1 and x2 y2 z2 should
%	be integer numbers, this program extends its limit to all
%	real numbers. If any of them are floating numbers, you
%	should specify how many digits of decimal that you would
%	like to preserve. Be aware that the length of output X Y
%	Z coordinates will increase in 10 times for each decimal
%	digit that you want to preserve. By default, the precision
%	is 0, which means that they will be rounded to the nearest
%	integer.
%
%  X	- a set of x coordinates on Bresenham's line
%
%  Y	- a set of y coordinates on Bresenham's line
%
%  Z	- a set of z coordinates on Bresenham's line
%
%  Therefore, all points in XYZ set (i.e. P(i) = [X(i) Y(i) Z(i)])
%  will constitute the Bresenham's line between P1 and P1.
%
%  Example:
%	P1 = [12 37 6];     P2 = [46 3 35];
%	[X Y Z] = bresenham_line3d(P1, P2);
%	figure; plot3(X,Y,Z,'s','markerface','b');
%
%  B.Pendleton.  line3d - 3D Bresenham's (a 3D line drawing algorithm)
%  ftp://ftp.isc.org/pub/usenet/comp.sources.unix/volume26/line3d, 1992
%
%  Also referenced by:
%
%  Fischer, J., A. del Rio (2004).  A Fast Method for Applying Rigid
%  Transformations to Volume Data, WSCG2004 Conference.
%  http://wscg.zcu.cz/wscg2004/Papers_2004_Short/M19.pdf
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
*/

/*
   if ~exist('precision','var') | isempty(precision) | round(precision) == 0
      precision = 0;
      P1 = round(P1);
      P2 = round(P2);
   else
      precision = round(precision);
      P1 = round(P1*(10^precision));
      P2 = round(P2*(10^precision));
   end

   d = max(abs(P2-P1)+1);
   X = zeros(1, d);
   Y = zeros(1, d);
   Z = zeros(1, d);

   x1 = P1(1);
   y1 = P1(2);
   z1 = P1(3);

   x2 = P2(1);
   y2 = P2(2);
   z2 = P2(3);

   dx = x2 - x1;
   dy = y2 - y1;
   dz = z2 - z1;

   ax = abs(dx)*2;
   ay = abs(dy)*2;
   az = abs(dz)*2;

   sx = sign(dx);
   sy = sign(dy);
   sz = sign(dz);

   x = x1;
   y = y1;
   z = z1;
   idx = 1;

   if(ax>=max(ay,az))			% x dominant
      yd = ay - ax/2;
      zd = az - ax/2;

      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;

         if(x == x2)		% end
            break;
         end

         if(yd >= 0)		% move along y
            y = y + sy;
            yd = yd - ax;
         end

         if(zd >= 0)		% move along z
            z = z + sz;
            zd = zd - ax;
         end

         x  = x  + sx;		% move along x
         yd = yd + ay;
         zd = zd + az;
      end

   else if (ay>=max(ax,az))		% y dominant
      xd = ax - ay/2;
      zd = az - ay/2;

      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;

         if(y == y2)		% end
            break;
         end

         if(xd >= 0)		% move along x
            x = x + sx;
            xd = xd - ay;
         end

         if(zd >= 0)		% move along z
            z = z + sz;
            zd = zd - ay;
         end

         y  = y  + sy;		% move along y
         xd = xd + ax;
         zd = zd + az;
      end

   else if (az>=max(ax,ay))		% z dominant
      xd = ax - az/2;
      yd = ay - az/2;

      while(1)
         X(idx) = x;
         Y(idx) = y;
         Z(idx) = z;
         idx = idx + 1;

         if(z == z2)		% end
            break;
         end

         if(xd >= 0)		% move along x
            x = x + sx;
            xd = xd - az;
         end

         if(yd >= 0)		% move along y
            y = y + sy;
            yd = yd - az;
         end

         z  = z  + sz;		% move along z
         xd = xd + ax;
         yd = yd + ay;
      end
   end

   if precision ~= 0
      X = X/(10^precision);
      Y = Y/(10^precision);
      Z = Z/(10^precision);
   end

   return;	
}
*/

// ****************************************************************************


// ***********************************************************************************************
/*
"A New Algorithm for Scan Conversion of a General Ellipse"
 ( 2002, C. Bond )
*/
/*
void
BondGeneralEllips( const C3Vector& in_coneAxis, const C3Vector& in_origin, const double& in_comptonAngle, 
	const double& in_sliceZ, int solveWhat, const double& knownCoordinateVal )
{
	double lambda = cos(in_comptonAngle);
	double lambda2 = lambda*lambda;
	double inverselambda2 = 1.0/lambda2;

	double nx = in_coneAxis.GetX();
	double ny = in_coneAxis.GetY();
	double nz = in_coneAxis.GetZ();
	double nx2 = nx * nx;
	double ny2 = ny * ny;
	double nz2 = nz * nz;

	double x1 = in_origin.GetX();
	double y1 = in_origin.GetY();
	double z1 = in_origin.GetZ();


//		double A = lambda2 - nx2;
//		double B = lambda2 - ny2;
//		double C = -1.0 * nx * ny;
//		double F = (lambda2 - nz2) * in_sliceZ;
		
//		A = A/(lambda2);
//		A = A*A;
//		B = B/(lambda2);
//		B = B*B;
//		C = in_sliceZ * C/(lambda2);
//		C = C*C;
//		F = in_sliceZ * F/(lambda2);
//		F = F*F;
//		
//		cout << "lambda: " << lambda << " lambda^2: " << lambda*lambda << " A:" << A << " B:" << B << " C:" << C << " F:" << F << endl;

	double zs = in_sliceZ;
	double zterm = (zs-z1);
	double zterm2 = zterm*zterm;

	// complete Wilderman
	// Factor for x^2
	double factor_x2 = 1.0 - 1.0 * inverselambda2 * (nx2);
	
	// Factor for y^2
	double factor_y2 = 1.0 - 1.0 * inverselambda2 * (ny2);

	// Factor for xy
	double factor_xy = -1.0 * inverselambda2 * (2.0 * nx * ny);

	// Factor for x
	double factor_x = -2.0 * x1;
	double inbrackets = ( nx2 * (-2.0 * x1) + 2.0 * nx * ny * (-1.0 * y1) + 2.0 * nx * nz * zterm );
	factor_x += -1.0 * inverselambda2 * inbrackets;

	// Factor for y
	double factor_y = -2.0 * y1;
	inbrackets = ( ny2 * (-2.0 * y1) + 2.0 * nx * ny * (-1.0 * x1) + 2.0 * ny * nz * zterm );
	factor_y += -1.0 * inverselambda2 * inbrackets;

	// f***ing rest
	double rest = x1*x1 + y1*y1 + zterm2;
	inbrackets = nx2*x1*x1 + ny2*y1*y1 + nz2*zterm2 + 2*nx*ny*x1*y1 + 2*nx*nz*(-1.0*x1*zs + x1*z1) + 2*ny*nz*(-1.0*y1*zs+y1*z1);
	rest += -1.0 * inverselambda2 * inbrackets;

	// cout
	cout << "factor_x2 = A: " << factor_x2 << "  ";
	cout << "factor_y2 = B: " << factor_y2 << "  ";
	cout << "factor_xy = C: " << factor_xy << "  ";
	cout << "factor_x = D: " << factor_x << "  ";
	cout << "factor_y = E: " << factor_y << "  ";
	cout << "rest = F: " << rest << endl;

	if (solveWhat == 0) // solve X, so Y,Z is known
	{
		// factor_x2 okay
		// factor_x = factor_x + factor_xy * knownCoordinateVal
		factor_x += factor_xy * knownCoordinateVal;
		// rest = rest + factor_y * knownCoordinateVal + factor_y2 * knownCoordinateVal^2
		rest += factor_y * knownCoordinateVal + factor_y2 * knownCoordinateVal*knownCoordinateVal;

		double xsol1, xsol2;
		solveQuadraticEquation( factor_x2, factor_x, rest, xsol1, xsol2);
		cout << "known y: " << knownCoordinateVal << ", x1+xsol1: " << x1+xsol1 << ", x1+xsol2: " << x1+xsol2 << endl;
	}
	else if (solveWhat == 1) // solve Y, so X,Z is known
	{
		// factor_y2 okay
		// factor_y = factor_y + factor_xy * knownCoordinateVal
		factor_y += factor_xy * knownCoordinateVal;
		// rest = rest + factor_x * knownCoordinateVal + factor_x2 * knownCoordinateVal^2
		rest += factor_x * knownCoordinateVal + factor_x2 * knownCoordinateVal*knownCoordinateVal;

		double ysol1, ysol2;
		solveQuadraticEquation( factor_y2, factor_y, rest, ysol1, ysol2);
		cout << "known x: " << knownCoordinateVal << ", y1+ysol1: " << y1+ysol1 << ", y1+ysol2: " << y1+ysol2 << endl;
	}
	else if (solveWhat == 2) // solve Z, so X,Y is known
	{
		// not possible yet...
		exit(1);
	}
}
*/


// ***********************************************************************************************

/*
"A New Algorithm for Scan Conversion of a General Ellipse"
 ( 2002, C. Bond )
*/

void
BondGeneralEllips( const C3Vector& in_coneAxis, const C3Vector& in_origin, const double& in_comptonAngle, 
	const double& in_sliceZ, int solveWhat, const double& knownCoordinateVal )
{
	double lambda = cos(in_comptonAngle);
	double lambda2 = lambda*lambda;
	double inverselambda2 = 1.0/lambda2;

	double nx = in_coneAxis.GetX();
	double ny = in_coneAxis.GetY();
	double nz = in_coneAxis.GetZ();
	double nx2 = nx * nx;
	double ny2 = ny * ny;
	double nz2 = nz * nz;

	double x1 = in_origin.GetX();
	double y1 = in_origin.GetY();
	double z1 = in_origin.GetZ();


//		double A = lambda2 - nx2;
//		double B = lambda2 - ny2;
//		double C = -1.0 * nx * ny;
//		double F = (lambda2 - nz2) * in_sliceZ;
		
//		A = A/(lambda2);
//		A = A*A;
//		B = B/(lambda2);
//		B = B*B;
//		C = in_sliceZ * C/(lambda2);
//		C = C*C;
//		F = in_sliceZ * F/(lambda2);
//		F = F*F;
//		
//		cout << "lambda: " << lambda << " lambda^2: " << lambda*lambda << " A:" << A << " B:" << B << " C:" << C << " F:" << F << endl;

	double zs = in_sliceZ;
	double zterm = (zs-z1);
	double zterm2 = zterm*zterm;

	// complete Wilderman
	// Factor for x^2
	double factor_x2 = 1.0 - 1.0 * inverselambda2 * (nx2);
	
	// Factor for y^2
	double factor_y2 = 1.0 - 1.0 * inverselambda2 * (ny2);

	// Factor for xy
	double factor_xy = -1.0 * inverselambda2 * (2.0 * nx * ny);

	// Factor for x
	double factor_x = -2.0 * x1;
	double inbrackets = ( nx2 * (-2.0 * x1) + 2.0 * nx * ny * (-1.0 * y1) + 2.0 * nx * nz * zterm );
	factor_x += -1.0 * inverselambda2 * inbrackets;

	// Factor for y
	double factor_y = -2.0 * y1;
	inbrackets = ( ny2 * (-2.0 * y1) + 2.0 * nx * ny * (-1.0 * x1) + 2.0 * ny * nz * zterm );
	factor_y += -1.0 * inverselambda2 * inbrackets;

	// f***ing rest
	double rest = x1*x1 + y1*y1 + zterm2;
	inbrackets = nx2*x1*x1 + ny2*y1*y1 + nz2*zterm2 + 2*nx*ny*x1*y1 + 2*nx*nz*(-1.0*x1*zs + x1*z1) + 2*ny*nz*(-1.0*y1*zs+y1*z1);
	rest += -1.0 * inverselambda2 * inbrackets;

	// cout
	cout << "factor_x2 = A: " << factor_x2 << "  ";
	cout << "factor_y2 = B: " << factor_y2 << "  ";
	cout << "factor_xy = C: " << factor_xy << "  ";
	cout << "factor_x = D: " << factor_x << "  ";
	cout << "factor_y = E: " << factor_y << "  ";
	cout << "rest = F: " << rest << endl;

	if (solveWhat == 0) // solve X, so Y,Z is known
	{
		// factor_x2 okay
		// factor_x = factor_x + factor_xy * knownCoordinateVal
		factor_x += factor_xy * knownCoordinateVal;
		// rest = rest + factor_y * knownCoordinateVal + factor_y2 * knownCoordinateVal^2
		rest += factor_y * knownCoordinateVal + factor_y2 * knownCoordinateVal*knownCoordinateVal;

		double xsol1, xsol2;
		solveQuadraticEquation( factor_x2, factor_x, rest, xsol1, xsol2);
		cout << "known y: " << knownCoordinateVal << ", x1+xsol1: " << x1+xsol1 << ", x1+xsol2: " << x1+xsol2 << endl;
	}
	else if (solveWhat == 1) // solve Y, so X,Z is known
	{
		// factor_y2 okay
		// factor_y = factor_y + factor_xy * knownCoordinateVal
		factor_y += factor_xy * knownCoordinateVal;
		// rest = rest + factor_x * knownCoordinateVal + factor_x2 * knownCoordinateVal^2
		rest += factor_x * knownCoordinateVal + factor_x2 * knownCoordinateVal*knownCoordinateVal;

		double ysol1, ysol2;
		solveQuadraticEquation( factor_y2, factor_y, rest, ysol1, ysol2);
		cout << "known x: " << knownCoordinateVal << ", y1+ysol1: " << y1+ysol1 << ", y1+ysol2: " << y1+ysol2 << endl;
	}
	else if (solveWhat == 2) // solve Z, so X,Y is known
	{
		// not possible yet...
		exit(1);
	}
}
