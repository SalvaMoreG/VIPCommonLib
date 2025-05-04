
// https://mathworld.wolfram.com/Circle-LineIntersection.html

#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{
	double x1, y1, z1(0.0);
	double x2, y2, z2(0.0);

	x1 = -120.0; y1 = 0.0;
	x2 =  120.0; y2 = 0.0;

	double radius(50.0);

	ofstream outfile("test_COINCS.dat_NEW");

	for (int ievent = 0; ievent < 10; ievent++)
	{
		double dx = x2 - x1;
		double dy = y2 - y1;
		double dr = sqrt(dx*dx + dy*dy);
	
		double Det = x1 * y2 - x2 * y1; // determinant
	
		double Discriminant = radius*radius * dr*dr - Det*Det;
	
		double distance = -1;
		if (Discriminant <= 0)
		{
			cout << "no solution found" << endl;
		}
		else
		{
			int sign = (dy >= 0) ? 1 : -1;
			double pm_x = sign * dx * sqrt(Discriminant);
			double o_x1 = (Det * dy + pm_x)/(dr*dr);
			double o_x2 = (Det * dy - pm_x)/(dr*dr);
	
			double pm_y = fabs(dy) * sqrt(Discriminant);
			double o_y1 = (Det * dx + pm_y)/(dr*dr);
			double o_y2 = (Det * dx - pm_y)/(dr*dr);
	
			double doy = o_y2 - o_y1;
			double dox = o_x2 - o_x1;
			distance = sqrt(doy*doy + dox*dox);
			cout << "A: " << o_x1 << " , " << o_y1 
			     << " B: " << o_x2 << " , " << o_y2 
				 << " distance: " << distance << endl;
		}

		outfile << z1 << " " << y1 << " " << x1 << "   " << 511 << "     "
		     << z2 << " " << y2 << " " << x2 << "   " << 511 << "     "
			 << distance << endl;

		y1 = y1 + 5;
		y2 = y2 + 5;
	}	
}
	
	
