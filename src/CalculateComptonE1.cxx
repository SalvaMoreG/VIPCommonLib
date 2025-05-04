
#include "CVIPUtils.h"
#include "CVIP3Vector.h"
#include "CVIPBaseCone.h"

#include "CVIPWildermanUtils.h"

#include <iostream>

double m_E_source( mass_electron_keV );

using namespace std;

int main()
{
	double angleDeg(1.0);
	double E1(0.0);
	CVIPBaseCone aCone;	
	C3Vector srcposition(0.0, 0.0, 0.0);
	double dx, dy, z(0.0);	
	C3Vector hitposition1(10.0 , 0.0 , z );
	C3Vector hitposition2; // (10.0 , 0.0 , 0.0 );	
	double length = 10.0; // hitposition2.GetLength();
	double sliceZ( z );
	
	std::vector<unsigned int> tmp_fovidxVec;
	CVIPFieldOfView fieldOfView;
	//	fieldOfView.SetFieldOfView( 61, -10.25, 20.25, 41, -10.25, 10.25, 5, -1.25, 1.25 );
	fieldOfView.SetFieldOfView( 21, -5.5, 15.5, 11, -5.5, 5.5, 1, -0.5, 0.5 );	
	// 	fieldOfView.SetFieldOfView( 10, -5.0, 15.0, 5, -5.0, 5.0, 1, -0.5, 0.5 );	
	C3Vector coords, minbounds, maxbounds;

	cout << "Total energy of gamma (< 0 = 511)" << endl;
	double tmpE;
	cin >> tmpE;
	if (tmpE > 0 )
		m_E_source = tmpE;

	while (angleDeg > 0)
	{
		cout << "Give angle (degrees)" << endl;
		cin >> angleDeg;
		if (angleDeg > 0)
		{
			double angleRad = angleDeg * kPI/180.0;
			bool ok = VIPUtils::GetComptonEnergyFromAngleRad(angleRad, m_E_source, E1 );
			cout << "E1: " << E1 << endl;

			// calculate hit position 2....
			dx = length * cos(angleRad);
			dy = length * sin(angleRad);
			cout << "dx: " << dx << " dy: " << dy << endl;
			hitposition2 = hitposition1 + C3Vector( dx, dy, 0 );
			cout << "HIT position 1: " << hitposition1 << endl;
			cout << "HIT position 2: " << hitposition2 << endl;			
			
			int error = aCone.Set( hitposition1, E1, hitposition2, (m_E_source - E1), m_E_source );
			if (!error)
			{
				cout << "CONE: " << endl;
				cout << "Angle: " << aCone.GetComptonAngle() << " rad; " 
								<< aCone.GetComptonAngle()*180.0/kPI << " degrees" << endl;
				cout << "Axis: " << aCone.GetComptonAxisDirection() << endl;
				cout << "Origin: " << aCone.GetComptonAxisOrigin() << endl;
			}
			else
			{
				cout << "ERROR! Invalid energy for Compton angle" << endl;
			}
			
			int scenario = WildermanMarch( sliceZ
						, aCone.GetComptonAxisOrigin(), aCone.GetComptonAxisDirection(), aCone.GetComptonAngle()
						, fieldOfView, tmp_fovidxVec, true );

			cout << "SCENARIO: " << scenario << endl;
			cout << "FOV bins on cone: " << endl;
			for (int i = 0; i < tmp_fovidxVec.size(); i++)
			{
				int index = tmp_fovidxVec[i];
				fieldOfView.GetVoxelCentre( index, coords);
				fieldOfView.GetVoxelBoundaries( index, minbounds, maxbounds);
				// cout << "[" << i << "], [" << index << "]: " << coords << endl;
				cout << "[" << i << "], [" << index << "]: " << coords 
                     << " between: [" << minbounds << " ; " << maxbounds << "]" << endl;
			}
			ShowFovVector( fieldOfView, tmp_fovidxVec );
		}
	}
}


