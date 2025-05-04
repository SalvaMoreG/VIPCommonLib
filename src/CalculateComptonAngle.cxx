
#include "CVIPUtils.h"
#include "CVIP3Vector.h"
#include "CVIPBaseCone.h"

#include <iostream>

using namespace std;

double m_E_source( mass_electron_keV );

int main()
{
	double E1(1.0);
	CVIPBaseCone aCone;
	C3Vector position1( 0.0, 0.0, 0.0);
	C3Vector position2(10.0, 0.0, 0.0);
	
	cout << "Total energy of gamma (< 0 = 511)" << endl;
	double tmpE;
	cin >> tmpE;
	if (tmpE > 0 )
		m_E_source = tmpE;

	while (E1 > 0)
	{
		cout << "Give Energy E1 (< 0 = stop)" << endl;
		cin >> E1;
		if (E1 > 0)
		{
			aCone.Set( position1, E1, position2, m_E_source - E1, m_E_source );
			double angle = aCone.GetComptonAngle();
			cout << "Angle: " << angle << " rad = " << angle*180.0/kPI << " degrees" << endl;

			bool isValid = VIPUtils::GetComptonAngleDegrees(E1, m_E_source, angle);
			if (!isValid) cout << "ERROR: invalid E1: " << E1 << endl;
			cout << "CHECK Angle: " << angle << " degrees" << endl;

		}
	}
}

