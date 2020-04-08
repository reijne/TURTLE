#include "../../../config.h"

using namespace std;

#include "PlainIndex.h"
#include "StoneyTwoElectronAccess.h"
#include "StoneyToPlainIndex.h"

#include <iostream>


int	main()
{
/*
StoneyTwoElectronAccess	test("/tmp/hanrath/Integrale/O3/fort.31", 10);

	test.putInList(PlainIndex(2));
	test.putInList(PlainIndex(5));
	test.putInList(PlainIndex(8));
	test.putInList(PlainIndex(1));
	test.putInList(PlainIndex(0));
	test.putInList(PlainIndex(6));
	
StoneyTwoElectronType	*p = 
		new StoneyTwoElectronType[test.getNumberOfListEntries()];
	test.getRequestedIntegrals(p);


	for ( INT i=0 ; i<test.getNumberOfListEntries() ; i++ )
		cout << i << ". " << p[i] << endl;

	
	delete p;
*/	

StoneyToPlainIndex	test("/tmp/hanrath/Integrale/O3/fort.31");

TwoElectronIntegralIndex<MOType>	ind(2, 2, 2, 1);
	test.set(ind);

/*
INT	i, j, k, l;
	while ( scanf("%d %d %d %d", &i, &j, &k, &l) )
	{	ind = TwoElectronIntegralIndex<MOType>(i+1, j+1, k+1, l+1);
		printf("(%d %d|%d %d)", i+1, j+1, k+1, l+1);
		test.set(ind);
	}
*/
	cout << test;
/*
INT	n;
	for ( INT i=0 ; i<10 ; i++ )
	{	n = iover2(i);
		cout << n << endl;
	}
*/
}
