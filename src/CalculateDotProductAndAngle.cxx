
#include "CVIP3Vector.h"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
	bool debug(false);
	std::string tmpStr;
	if (argc > 1)
	{
		tmpStr = argv[1];
		if (tmpStr == "-debug")
		{
			debug = true;
		}
	}
	cout << "Give 1st vector coordinates: " << endl;
	C3Vector v1; cin >> v1;
	cout << "1st vector: " << v1 << endl;

	cout << "Give 2nd vector coordinates: " << endl;
	C3Vector v2; cin >> v2;
	cout << "2nd vector: " << v2 << endl;

	double angleRad = v1.GetScalarProductAngleRadians(v2);
	double dotproduct = v1*v2;

	cout << "v1*v2: " << v1*v2 
         << " angle: " << angleRad << " radians = " 
         << angleRad/kPI << "pi = " 
         << angleRad*180.0/kPI << " degrees" 
		 << " = " << 180.0 - angleRad*180.0/kPI << " degrees" 
         << endl;

	if (debug)
	{
		cout << " |v1|: " << v1.GetLength()
			 << " |v2|: " << v2.GetLength()
			 << " |v1|*|v2|: " << v1.GetLength() * v2.GetLength()
	   	  	 << " v1*v2/(|v1|*|v2|): " << v1*v2/( v1.GetLength() * v2.GetLength())
			 << " acos...: " << acos(v1*v2/( v1.GetLength() * v2.GetLength()))
			 << endl;
	}

	return 0;
}


