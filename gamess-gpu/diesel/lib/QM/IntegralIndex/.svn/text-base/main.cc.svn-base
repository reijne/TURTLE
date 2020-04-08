#include "../../../config.h"

using namespace std;

#include "../MRTree/Container/IndexedContainer.h"

#include "TwoElectronIntegralIndex.h"
#include "TwoElectronIntegralTriadeIndex.h"
#include "TwoElectronIntegralIndexTemplate.h"
#include "TwoElectronIntegralCbExIndex.h"
#include "IndexMask.h"
#include "../MO/MRMOs.h"

#include <iostream>


int main()
{
MRMOs	mrmos("/tmp/hanrath/Integrale/O3/fort.31");
	GeneralizedMO::mrmos = &mrmos;
	
	
IndexedContainer<double>	integrals(1000);
/*
	integrals[3] = 4;
	cout << integrals[3] << endl;


*/


/*
TwoElectronIntegralIndex<MOType> t(5, 2, 4, 6);
//TwoElectronIntegralIndex<MOType> t(10, 5, 1, 3);
/TwoElectronIntegralIndex<MOType> t(5, 6, 9, 8);
	cout << t << endl;

TwoElectronIntegralIndex<MOType> t1(5, 4, 2, 6);
//TwoElectronIntegralIndex<MOType> t1(10, 1, 5, 3);
//TwoElectronIntegralIndex<MOType> t1(5, 8, 6, 9);
	cout << t1 << endl;

TwoElectronIntegralTriadeIndex tt(t);
	cout << tt << endl;
	tt = TwoElectronIntegralTriadeIndex(t1);
	cout << tt << endl;

TwoElectronIntegralCbExIndex tt1(t, t1);
	cout << tt1 << endl;

//	exit(0);

//	for ( INT i=0 ; i<20000000 ; i++ )
//		tt.set(t);
//		tt = TwoElectronIntegralTriadeIndex(t);
	for ( INT i=0 ; i<20000000 ; i++ )
//		tt1 = TwoElectronIntegralCbExIndex(t, t1);
		tt1.set(t, t1);
	exit(0);
*/

TwoElectronIntegralIndex<GeneralizedMO> tgen(2, 1, 3, 4);
	tgen[0].setType(GeneralizedMO::IrRep_Space);
	tgen[0].setSignature(0);
	tgen[1].setType(GeneralizedMO::IrRep_Space);
	tgen[1].setSignature(1);
	cout << "tgen= " << tgen << endl;

MOType	MOs[2] = { 23, 88 };

TwoElectronIntegralIndexTemplate temp(1, 2, 3, 4);
	temp = tgen;
	cout << "temp= " << temp << endl;

//	temp.setMaskK(0);
	for ( INT i=0 ; i<100000000 ; i++ )
		temp.setFromTemplate(MOs);

	cout << temp << endl;

/*
IndexMask	im(4);
	im.setMask(t, 2);
	cout << im << endl;
	
	t.setFromMask(im, 9);
	cout << t << endl;
*/
}
