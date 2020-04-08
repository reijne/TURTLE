#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

#include "../../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../../lib/QM/MO/MOMapping.h"
#include "../../../../lib/QM/MO/MOIrReps.h"
#include "../../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../../lib/Container/DiskBuffer.h"
#include "../../../../lib/QM/RepresentationMatrices/CIVectors.h"


#include "../../../../config.h"
#include "../../../../VersionDate.h"

int	main(INT argc, char **argv)
{
	cout << "lstconfs (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc<2 || argc>4 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program prints the configurations in a tree (\"ConfTree.dat\"-format)." << endl; 
		cerr << endl;
		cerr << "usage: lstconfs \"ConfTree.dat\"-file [\"Eigenvectors.dat\"-file] [\"fort.31\"-file]" << endl;
		cerr << endl;
		exit(1);
	}
	
ifstream	ConfIn(argv[1]);
MOMapping	dummyMoMapping(ConfIn);


typedef double VectorType;

Fort31File	f31(argc == 4 ? argv[3] : "fort.31");
MOIrReps	moirreps(f31);

NExternalsDiag	MRDCIConfs(ConfIn, &moirreps);

DiskBuffer	*evBuf = NULL;
CIVectors<VectorType>	*ev = NULL;
INT	SAFs = MRDCIConfs.getNumberOfTotalSpinAdaptedFunctions();
INT	nRoots = 0;
	if ( argc>=3 )
	{
		evBuf = new DiskBuffer(argv[2]);
		nRoots = evBuf->getNumberOfObjects();
		ev = new CIVectors<VectorType>(SAFs, nRoots);

		for ( INT i=0 ; i<evBuf->getNumberOfObjects() ; i++ )
			evBuf->get(i, ev->getP(i));
		delete evBuf;
	}

INT	ii = 0;
	for ( MRTreeIterator i = MRDCIConfs.firstInTree() ;
			!MRDCIConfs.isLastInTree(i) ; MRDCIConfs.nextInTree(i))
	{
		cout << MRDCIConfs.getConfigurationSAFNr(i);
		if ( ev )
			for ( INT j=0 ; j<nRoots ; j++ )
				cout << "\t" << (*ev)(ii, j);
		ii++;
		cout << endl;
	}
	
	if ( argc==3 )
		delete ev;
	
	return 0;
}
