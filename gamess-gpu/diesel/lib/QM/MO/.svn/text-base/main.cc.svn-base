#include "../../../config.h"

using namespace std;

#include "GeneralizedMO.h"

#include "MRMOs.h"

#include <iostream>


int	main()
{
/*GeneralizedMO	mo1(2, 3);



	exit(0);
*/	
	
MRMOs	mrmos("/tmp/hanrath/Integrale/O2/fort.31");


	for ( INT i=0 ; i<12 ; i++ )
		mrmos.setInternal(i);

	mrmos.initIntExt();
	cout << mrmos << endl;
	
	
	
	GeneralizedMO::mrmos = &mrmos;
GeneralizedMO	mo2(2), mo1;
//	mo1.set(20);
//	mo2.set(21);
	mo1.set(19, GeneralizedMO::IrRep_Space);
	cout << mo1 << endl;
	cout << mo2 << endl;

/*	mo2 = mo1;

	cout << mo1 << endl;
	cout << mo2 << endl;
*/
	cout << (mo1==mo2) << endl;
	cout << (mo1<mo2) << endl;
	cout << (mo1>mo2) << endl;
}
