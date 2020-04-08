#include "../../../config.h"

#include "FourIndexIntegralContainer.h"
#include "ExtPosContainer.h"
#include "NTupelContainer.h"

#include "../IO/Fortran/Fort31File.h"


int	main()
{
/*INT	sub[4];
NTupelContainer	t(1000, sub);

TwoElectronIntegralTriadeIndex	ind(3, 4, 5, 6, 1);
	printf("%f", t[ind]);
*/	

Fort31File	f31;
	f31.name = "/tmp/hanrath/Integrale/C6/fort.31";
	f31.format = Fort31RecordFormatNew;

MRMOs	mrmos(f31);
	mrmos.initIntExt();
FourIndexIntegralContainer	cont(mrmos);


	cout << "total number of containers: " << cont.getNumberOfLeaves() << endl;
	cout << "total number of integrals: " << cont.getNumberOfIntegrals() << endl;

	cout << "allocating memory for integrals..." << endl;
	
	if ( cont.allocate() )
		cout << "success" << endl;
	else
		cout << "failed" << endl;
	


	cont.loadIntegrals(f31);
	
	getchar();



/*
TwoElectronIntegralTriadeIndex	t(3, 4, 5, 6, 1);

	printf("%f\n", cont[t]);
	
	
double	v = 76, h = 8;
*/

/*
SymmetryContainerIJKL	test(100000);
TwoElectronIntegralTriadeIndex	t(3, 4, 5, 6, 1);
INT	ind[4] = {1, 4, 3, 1};
INT	mult[4] = {8, 3, 2, 1};


//double	v = 76, h = 8;
*/
/*
	printf("%d\n", (ind[0]*mult[0] +
			ind[1]*mult[1] +
			ind[2]*mult[2] +
			ind[3]));
*/


/*	for ( INT i=0 ; i<10000000 ; i++ )
	{	t[0] = i & 1;
		t[1] = i & 1;
		t[2] = i & 15;
		v = cont[t];
		h *= v;
	}
*/

/*	for ( INT i=0 ; i<10000000 ; i++ )
	{	ind[0] = i & 1;
		ind[1] = i & 1;
		ind[2] = i & 15;
		test.get(ind, 2, v);
		h *= v;
	}
*/
//	printf("%f\n", h);

}
