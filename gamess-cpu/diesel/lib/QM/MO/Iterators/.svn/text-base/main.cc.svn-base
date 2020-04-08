#include "../../../config.h"

#include "MORestriction.h"
#include "MORestrictionState.h"
#include "../../Configuration/Configuration.h"

#include "MOEquivalenceProjected.h"


/*#include "OrderMatrixIterator.h"

#include "TupelIterator.h"

#include "MatchingMOListIterator.h"
*/


#include "MOIterator.h"
#include "MOAccess.h"


int	main()
{
/*
OrderMatrixIterator	iter(1, 1);


	iter.first();
	while ( !iter.isEnd() )
	{
		cout << iter;
		iter.next();
	}

*/
/*
TupelIterator	tupel(8);

	tupel.first();
	while ( !tupel.isEnd() )
	{
		cout << tupel << endl;
		tupel.next();
	}

*/
/*
Fort31File f31("/tmp1/hanrath/calc/SiCluster/Si6H6/BenzvalenDoppel/1A1_avg/fort.31", Fort31RecordFormatNew);
MRMOs	*mrmos = new MRMOs(f31);

	cout << *mrmos << endl;
*/
//MatchingMOListIterator	match();


/*
MORestriction	moRestriction("1-10=2   12-14=1 20-30>=0 18,19<=1 40-50=0");
//MORestriction	moRestriction("1-10=1  40-50=0");

	moRestriction.setMaxMO(100);

	cout << moRestriction << endl;
	if ( moRestriction.isIllegal() )
		cout << "illegal!" << endl;
		
	cout << moRestriction.getRestriction() << endl;

Configuration<MOType>	conf;
	conf.create(2);
//	conf.create(2);
	conf.create(13);
	conf.create(18);
	conf.create(25);
	conf.create(26);
//	conf.create(44);

MORestrictionState	state(moRestriction, conf);

	cout << "impossible: " << state.isImpossible() << endl;
	cout << "check: " << state.check() << endl;

Configuration<MOType>	conf2;
	conf2.create(3);

MORestrictionState	state2(moRestriction, conf2);

	state |= state2;
	cout << "impossible: " << state.isImpossible() << endl;
	cout << "check: " << state.check() << endl;
*/

/*
MOType	p[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
MOEquivalenceProjected::distribute(10, p, 5);
*/


Fort31File f31("/tmp1/hanrath/calc/H2O/fort.31", Fort31RecordFormatNew);
MRMOs	*mrmos = new MRMOs(f31);

	mrmos->setInternal(1);
	mrmos->setInternal(2);
	mrmos->setInternal(3);
	mrmos->setInternal(4);
	mrmos->setInternal(9);
	mrmos->setInternal(10);
	mrmos->initIntExt();

	cout << *mrmos << endl << endl << endl;

/*
MOIterator	iter(2, mrmos, 1);

Configuration<MOType>	conf(iter);

	cout << "n=" << iter.getN() << endl;

	while ( !iter.isEnd() )
	{
		cout << iter << endl;
		iter.next();
	}
*/	
MOAccess	moAccess(2, mrmos);
	cout << moAccess(14, 13) << endl;
	cout << moAccess(13, 14) << endl;

}
