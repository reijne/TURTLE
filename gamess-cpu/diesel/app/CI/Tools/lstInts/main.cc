#include <iostream>
#include <stdio.h>

#include "../../../../lib/QM/Parallel/SharedMemory.h"
#include "../../../../lib/QM/MO/MRMOs.h"
#include "../../../../lib/QM/IntegralContainer/FourIndexIntegralContainer.h"
#include "../../../../lib/QM/IntegralContainer/TwoIndexIntegralContainer.h"
#include "../../../../lib/QM/IntegralIndex/OneElectronIntegralIndex.h"
#include "../../../../lib/QM/IntegralIndex/TwoElectronIntegralIndex.h"

#include "../../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../../lib/QM/IO/Fortran/Fort31FirstRecord.h"

#include "../../../../config.h"
#include "../../../../VersionDate.h"

using namespace std;

int	main(INT argc, char **argv)
{
	cout << "lstints (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc!=1 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program prints integrals from fort.31-file." << endl; 
		cerr << endl;
		exit(1);
	}
	
Fort31File	Fort31("fort.31");
MRMOs	mrmos(Fort31);
	mrmos.initIntExt();


FourIndexIntegralContainer	*intContainer = new FourIndexIntegralContainer(mrmos);

	INT	size = intContainer->getDataSize();

SharedMemory	*Int4sharedMem = new SharedMemory(size, SharedMemory::StandAlone);


		intContainer->setSharedMem(Int4sharedMem);


		cout << "total number of two electron integral containers: " << 
			intContainer->getNumberOfLeaves() << endl;
		cout << "total number of two electron integrals: " << 
			intContainer->getNumberOfIntegrals() << endl;

		cout << "allocating memory for two electron integrals..." << flush;

		if ( intContainer->allocate() )
				cout << "success" << endl;
		else
		{
			cout << "failed" << endl;
			exit(1);
		}

		intContainer->loadIntegrals(Fort31);


TwoIndexIntegralContainer	*int2Container = new TwoIndexIntegralContainer(mrmos);

	double	core = int2Container->loadIntegrals(Fort31);


//--------------------------------------------------------------------------

	Fort31FirstRecord	f31(Fort31);


		cout << "VNUC = " << f31.getVNUC() << endl;
		cout << "ZERO = " << f31.getZERO() << endl;
		cout << "CORE = " << core << endl;



OneElectronIntegralIndex	ind1;
TwoElectronIntegralIndex<MOType>	ind2;

INT	mo[4];
/*



double	E = 0;
const INT	nMOs = 5;
INT	MO[nMOs] = { 1, 2, 3, 11, 16 };


	for ( INT i=0 ; i<nMOs ; i++ )
	{
		ind1.set(MO[i], MO[i]);
		E += 2*(*int2Container)[ind1];
		cout << "i=" << i << ", E=" << E << endl;
		for ( INT j=0 ; j<nMOs ; j++ )
		{
			ind2.set(MO[i], MO[i], MO[j], MO[j]);
			E += 2*(*intContainer)[ind2];

			ind2.set(MO[i], MO[j], MO[i], MO[j]);
			E += (*intContainer)[ind2];
		}
		cout << "i=" << i << ", E=" << E << endl;
	}


	cout << "E=" << E << endl;
*/

	while ( cin )
	{
	char	buf[1000];
		cin.getline(buf, 1000);
	INT	n = sscanf(buf, "%d %d %d %d\r", &mo[0], &mo[1], &mo[2], &mo[3]);
		if ( n==2 )
		{
			ind1.set(mo[0], mo[1]);
			cout << "(" << mo[0] << " " << mo[1] << ") = " << (*int2Container)[ind1] << endl;
		}
		else
		if ( n==4 )
		{
		ind2.set(mo[0], mo[1], mo[2], mo[3]);
		cout << "(" << mo[0] << " " << mo[1] << "|" << mo[2] << " " << mo[3] << ") = " << (*intContainer)[ind2] << endl;
		}
		else
			break;
	}



}
