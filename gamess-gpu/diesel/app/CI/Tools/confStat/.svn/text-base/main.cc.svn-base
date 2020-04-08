#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;

#include "../../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../../lib/QM/MO/MOMapping.h"
#include "../../../../lib/QM/MO/MOIrReps.h"
#include "../../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"

#include "../../../../lib/Math/MatrixVector/Vector.h"
#include "../../../../lib/Container/DiskBuffer.h"

#include "../../../../lib/Container/GrepAwk.h"


#include "../../Selector/MRConfInput.h"

#include "../../../../lib/QM/MO/Iterators/MOEquivalence.h"

#include "../../../../lib/QM/Configuration/ConfigurationSAFNr.h"


#include "../../../../config.h"
#include "../../../../VersionDate.h"

int	main(INT argc, char **argv)
{
	cout << "confStat (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE 
		<< endl << endl;

	if ( argc<2 || argc>4 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program counts excitations relative to a certain configuration." << endl; 
		cerr << endl;
		cerr << "usage: confStat Threshold <\"selector input\" [\"fort.31\"-file]" << endl;
		cerr << endl;
		exit(1);
	}


String	Threshold(argv[1]);
ifstream	ConfIn("ConfTree.dat."+Threshold);
MOMapping	dummyMoMapping(ConfIn);


MRConfInput	mrinp;
	cin >> mrinp;

ConfigurationSet	refs;

	Pix j = mrinp.getRefConfSet().first();
	while ( j )
	{
		refs.add(mrinp.getRefConfSet()(j));
		mrinp.getRefConfSet().next(j);
	}



Fort31File	f31(argc == 4 ? argv[3] : "fort.31");
MOIrReps	moirreps(f31);

	cout << "reading tree..." << endl;
NExternalsDiag	tree(ConfIn, &moirreps);
	cout << "OK" << endl;

const INT nRoots = tree.getNumberOfRoots();


typedef double VectorType;

const INT maxRefs = refs.length();
const INT maxExc = 100;

INT	exc[maxRefs][maxExc];
	memset(exc, 0 , maxRefs*maxExc*sizeof(INT));
double	ci2[maxRefs][nRoots];
	memset(ci2, 0 , maxRefs*nRoots*sizeof(double));


	cout << "reading diag.out." << Threshold << "..." << flush;

ifstream	diagOut("diag.out."+Threshold);
GrepAwk	ga(diagOut);

	if ( !ga.grep("R e s u l t s") )
	{
		cout << "error: output contains no wave function." << endl;
		exit(1);
	}
	
	for ( INT r=0 ; r<nRoots ; r++ )
	{
		if ( !ga.grep("root # ") )
		{
			cout << "error: output contains no wave function." << endl;
			exit(1);
		}
		ga += 3;
		for ( INT i=0 ; i<refs.length() ; i++ )
		{
			while ( !ga.illegal() && ga.getWord(3)!="ref." )
				ga++;
			if ( ga.illegal() )
			{
				cout << "error: output contains no wave function." << endl;
				exit(1);
			}
			ci2[i][r] = atof(ga.getWord(1));
			ga++;
		}
		
		
	}
	

	cout << "OK" << endl;

	cout << "counting exciation orders..." << flush;


	for ( MRTreeIterator i = tree.firstInTree() ;
			!tree.isLastInTree(i) ; tree.nextInTree(i))
	{
	ConfigurationSAFNr<MOType> conf = tree.getConfigurationSAFNr(i);
		
	Pix	refPix = refs.first();
	INT	refi = 0;
		while ( refPix )
		{
			exc[refi][Configuration<MOType>::calcExcitationOrder(conf, refs(refPix))]++;
			refi++;
			refs.next(refPix);
		}
	}
	cout << "OK" << endl;



//	use eigenvectors (optional)
/*	cout << "reading eigenvectors..." << flush;
const INT SAFs = tree.getNumberOfTotalSpinAdaptedFunctions();

String	evName("Eigenvectors.dat."+Threshold);
ifstream	
Vector<VectorType>	*ev = new Vector<VectorType> [nRoots] (SAFs);
DiskBuffer	*evBuf = new DiskBuffer(SAFs*sizeof(VectorType),
		evName.chars(), 0);

	for ( INT i=0 ; i<nRoots ; i++ )
		evBuf->get(i, ev[i].getP());

	cout << "OK" << endl;

	cout << "counting exciation orders..." << flush;
const INT maxRefs = refs.length();
const INT maxExc = 100;

INT	exc[maxRefs][maxExc];
	memset(exc, 0 , maxRefs*maxExc*sizeof(INT));
double	ci2[maxRefs][nRoots];
		memset(ci2, 0 , maxRefs*nRoots*sizeof(double));

	for ( MRTreeIterator i = tree.firstInTree() ;
			!tree.isLastInTree(i) ; tree.nextInTree(i))
	{
	ConfigurationSAFNr<MOType> conf = tree.getConfigurationSAFNr(i);
		
	INT	refI = 0;
	Pix	refPix = refs.first();
	INT	refi = 0;
		while ( refPix )
		{
			exc[refi][Configuration<MOType>::calcExcitationOrder(conf, refs(refPix))]++;

			if ( conf.isReference() && conf==refs(refPix) )
				refI = refi;


			refi++;
			refs.next(refPix);
		}

		if ( conf.isReference() )
		{
			for ( INT r=0 ; r<nRoots ; r++ )
			{
			double	h = 0;
				for ( INT j=conf.getSAFNr() ; 
					j<conf.getSAFNr()+conf.getSAFInc() ; j++ )
					h += ev[r][j] * ev[r][j];
				ci2[refI][r] = h;
			}
		}
	}
	delete []  ev;
	delete evBuf;
	cout << "OK" << endl;
*/
INT i = 0;
	for ( i=maxExc-1 ; i>=0 ; i-- )
	{
	INT	breaked = 0;
		for ( INT j=0 ; j<maxRefs ; j++ )
			if ( (breaked=exc[j][i]) )
				break;
		if ( breaked )
			break;
	}


	cout << endl << endl << endl;

Pix	iRef = refs.first();
INT	k = 0;
	for ( INT j=1 ; j<=i ; j++ )
		cout << setw(10) << "order";
	cout << "   |   ";
	for ( INT j=0 ; j<nRoots ; j++ )
		cout << setw(6) << "root";
	cout << endl;
	for ( INT j=1 ; j<=i ; j++ )
		cout << setw(10) << j;
	cout << "   |   ";
	for ( INT j=0 ; j<nRoots ; j++ )
		cout << setw(6) << j+1;
	cout << "    |           configuration" << endl;
	for ( INT j=1 ; j<=i ; j++ )
		cout << "----------";
	cout << "---+---";
	for ( INT j=0 ; j<nRoots ; j++ )
		cout << "------";
	cout << "----+------------------------------------------------------------" << endl;
	cout.setf(ios::fixed);
	while ( iRef )
	{
		for ( INT j=1 ; j<=i ; j++ )
			cout << setw(10) << exc[k][j];
		cout << "   |   ";
		for ( INT j=0 ; j<nRoots ; j++ )
		{
			if ( ci2[k][j]>0.01 )
				cout << setw(6) << setprecision(2) << ci2[k][j];
			else
				cout << "  ----";
		}
		cout << "    |  #" << setw(4) << k+1 << ":   ";
		cout << refs(iRef) << endl;
		refs.next(iRef);
		k++;
	}
/*
INT i = 0;
	for ( i=maxExc-1 ; i>=0 ; i-- )
		if ( exc[i] )
			break;
	for ( INT j=0 ; j<=i ; j++ )
		cout << setw(10) << exc[j] << " of order " << j << endl;
*/		

	return 0;
}
