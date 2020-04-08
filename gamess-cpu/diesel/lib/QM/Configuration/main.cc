#include "../../../config.h"

#include "Configuration.h"
#include "TableCase.h"
#include "DiffConf.h"
#include "Excitation.h"
#include "../../Math/etc/BinomialCoefficient.h"
#include "../MO/GeneralizedMO.h"
#include "../MO/MRMOs.h"



//#define TMOTYPE GeneralizedMO
#define TMOTYPE MOType



int	main()
{
Configuration<MOType>	conf;

	cin >> conf;
	cout << conf;


/*
BinomialCoefficient	binom(MAXOPENSHELLS);

MRMOs	mrmos("/tmp/hanrath/Integrale/O2/fort.31");
	mrmos.initIntExt();
	cout << mrmos << endl;


GeneralizedMO::mrmos = &mrmos;

Configuration<TMOTYPE>	conf1(4, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6);
Configuration<TMOTYPE>	conf2(8, 4, 1, 2, 7, 8, 9, 10, 36, 37, 3, 4, 5, 6);
DiffConf<TMOTYPE>	dc;
TableCase<TMOTYPE>	tc(&binom);

	dc.calcDiffConf(conf1, conf2);
	cout << dc << endl;
	tc.calc(dc);
	cout << tc << endl;

*/
/*
Configuration<TMOTYPE>	conf1(4, 0, 7, 8, 9, 10);
Configuration<TMOTYPE>	conf2(4, 0, 8, 9, 10, 11);
DiffConf<TMOTYPE>	dc;

	dc.calcDiffConf(conf1, conf2);
	cout << dc << endl;

	conf1 = Configuration<TMOTYPE>(4, 0, 7, 8, 9, 10);
	conf2 = Configuration<TMOTYPE>(3, 0, 8, 9, 10);
	dc.calcDiffConf(conf1, conf2);
	dc.addExternal(
		Configuration<MOType>(),
		Configuration<MOType>(1, 0, 11));
	cout << dc << endl;
*/	

/*
Configuration<TMOTYPE>	conf = Configuration<TMOTYPE>(4, 5, 7, 8, 9, 10, 2, 3, 4, 5, 6);
DiffConf<TMOTYPE>	dc;
	dc.calcDiffConf(conf, conf);
	dc.addExternal(
		Configuration<MOType>(0, 1, 36),
		Configuration<MOType>(2, 0, 36, 37));
	cout << dc << endl;
*/

/*
Configuration<TMOTYPE>	conf1(2, 3, 3, 4, 1, 2, 5);
//	cout << conf1 << endl << endl;

Configuration<TMOTYPE>	conf2;
TableCase<TMOTYPE>	tc(&binom);
DiffConf<TMOTYPE>	dc;

	conf1 = Configuration<TMOTYPE>(2, 2, 1, 5, 3, 4);
	conf2 = Configuration<TMOTYPE>(4, 1, 1, 2, 3, 5, 4);
	dc.calcDiffConf(conf1, conf2);
	cout << dc << endl;
	conf2.excite(1, 6);
	dc.calcDiffConf(conf1, conf2);
	cout << dc << endl;
	cout << conf1 << endl;
	cout << conf2 << endl;
	tc.calc(conf1, conf2);
	cout << tc << endl;

*/
	exit(0);

Excitation	exc;

/*
	conf1 = Configuration<TMOTYPE>(2, 2, 7, 8, 1, 5);
	conf2 = Configuration<TMOTYPE>(2, 2, 2, 9, 1, 5);
	cout << conf1 << ", " << conf2 << endl;
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = Configuration<TMOTYPE>(2, 2, 6, 9, 1, 5);
	cout << conf1 << ", " << conf2 << endl;
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(4, 0, 1, 5, 7, 8);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 5, 9);
	cout << conf1 << ", " << conf2 << endl;
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = Configuration<TMOTYPE>(4, 0, 1, 5, 6, 9);
	cout << conf1 << ", " << conf2 << endl;
	tc.calc(conf1, conf2);
	cout << tc << endl;
	exit(0);
*/


/*
DiffConf<TMOTYPE>	dc;

//	conf1 = Configuration<TMOTYPE>(2, 2, 1, 4, 2, 3);
//	conf2 = Configuration<TMOTYPE>(4, 1, 2, 4, 6, 7, 5);
//	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 1);
//	conf2 = Configuration<TMOTYPE>(2, 1, 1, 2, 3);
//	conf1 = Configuration<TMOTYPE>(4, 0, 1, 3, 5, 15);
//	conf2 = Configuration<TMOTYPE>(4, 0, 1, 4, 5, 15);
	conf1 = Configuration<TMOTYPE>(4, 10, 7, 8, 34, 46, 1, 2, 3, 4, 5, 6, 9, 10, 35, 47);
	conf2 = Configuration<TMOTYPE>(3, 10, 8, 35, 46, 1, 2, 3, 4, 5, 6, 9, 10, 34, 47);
	cout << conf1 << endl;
	cout << conf2 << endl;
	cout << calcExcitationOrder(conf1, conf2);
//	exc.calcExcitation(conf1, conf2);
//	cout << exc << endl;
	dc.calcDiffConf(conf1, conf2);
//	cout << dc << endl;
	conf1 = Configuration<TMOTYPE>(0, 0);
	conf2 = Configuration<TMOTYPE>(1, 0, 11);
	dc.addExternal(conf1, conf2);
	cout << dc << endl;
	tc.calc(dc);
	cout << tc << endl;
	exit(0);
*/	


/*
DiffConf<GeneralizedMO>	dc;
TableCase<GeneralizedMO>	tcGen(&binom);

	mrmos.setInternal(1);
	mrmos.setInternal(2);
	mrmos.setInternal(4);

	conf1 = Configuration<MOType>(2, 0, 1, 4);
	conf2 = Configuration<MOType>(2, 0, 2, 4);
	cout << conf1 << endl;
	cout << conf2 << endl;
	dc.calcDiffConf(conf1, conf2);
	cout << dc << endl;
	conf1 = Configuration<MOType>(1, 0, 12);
	conf2 = Configuration<MOType>(1, 0, 11);
	dc.addExternal(conf1, conf2);
	cout << dc << endl;
	tcGen.calcLess3(dc);
	cout << tcGen << endl;
	exit(0);
*/
	

/*	conf1 = Configuration<TMOTYPE>(10, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1);
	conf2 = Configuration<TMOTYPE>(10, 1, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 3);
//	exc.calcExcitation(conf1, conf2);
	cout << exc << endl;
	for ( INT i=0 ; i<1000000 ; i++ )
		tc.calc(conf1, conf2);
//		exc.calcExcitation(conf1, conf2);
	cout << tc << endl;
	exit(0);
*/
	
	
/*	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 1);
	conf2 = Configuration<TMOTYPE>(2, 1, 5, 6, 7);
	tc.calc(conf1, conf2);
//	exc.calcExcitation(conf1, conf2);
	cout << exc << endl;
	cout << tc << endl;
	cout << isExcitationGreaterTwo(conf1, conf2) << endl;
	exit(0);
*/	


/*	
//	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 1);
//	conf2 = Configuration<TMOTYPE>(2, 1, 5, 6, 7);
	
	conf2 = conf1;
	conf2.excite(5, 6);
	conf2.excite(5, 6);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(5, 6);
	conf2.excite(5, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(5, 6);
	conf2.excite(5, 7);
	tc.calc(conf1, conf2);
	cout << tc << endl;
//	exit(0);

	conf2 = conf1;
	conf2.excite(5, 6);
	conf2.excite(2, 6);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(3, 6);
	conf2.excite(4, 7);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(5, 3);
	conf2.excite(2, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

//

	cout << "---------------------------------------------" << endl;

	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 1);
	conf2 = Configuration<TMOTYPE>(2, 1, 1, 2, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 2);
	conf2 = Configuration<TMOTYPE>(2, 1, 1, 2, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 1);
	conf2 = Configuration<TMOTYPE>(2, 1, 1, 2, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(2, 1, 3, 4, 2);
	conf2 = Configuration<TMOTYPE>(2, 1, 1, 2, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

//

	conf1 = Configuration<TMOTYPE>(1, 1, 4, 1);
	conf2 = Configuration<TMOTYPE>(3, 0, 1, 2, 3);
//	conf2.getOpenShellRef(1).setType(GeneralizedMO::IrRep_Space);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(1, 1, 4, 2);
	conf2 = Configuration<TMOTYPE>(3, 0, 1, 2, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(1, 1, 4, 3);
	conf2 = Configuration<TMOTYPE>(3, 0, 1, 2, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

//

	conf1 = Configuration<TMOTYPE>(1, 2, 4, 2, 3);
	conf2 = Configuration<TMOTYPE>(3, 1, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(1, 2, 4, 1, 3);
	conf2 = Configuration<TMOTYPE>(3, 1, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(1, 2, 4, 1, 2);
	conf2 = Configuration<TMOTYPE>(3, 1, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

//

	conf1 = Configuration<TMOTYPE>(0, 2, 1, 2);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(0, 2, 1, 3);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(0, 2, 1, 4);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(0, 2, 2, 3);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(0, 2, 2, 4);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf1 = Configuration<TMOTYPE>(0, 2, 3, 4);
	conf2 = Configuration<TMOTYPE>(4, 0, 1, 2, 3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;


	cout << "--------------------------------" << endl;

	conf1 = Configuration<TMOTYPE>(2, 3, 3, 4, 1, 2, 5);
//----------------------------------------------------

	conf2 = conf1;
	conf2.excite(3, 6);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(5, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(4, 3);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	conf2 = conf1;
	conf2.excite(3, 4);
	tc.calc(conf1, conf2);
	cout << tc << endl;

	cout << "--------------------------------" << endl;
//----------------------------------------------------

	conf2 = conf1;
	tc.calc(conf1, conf2);
	cout << tc << endl;
*/
}


