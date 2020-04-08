#include "NExternalsDiag.h"
#include "MRCIMatrix.h"


#include "MRTreeIterator.h"
#include "MOMapping.h"

#include "ConfigurationSet.h"

#include "StartUp.h"
#include "VarSizeReadOnlyCache.h"


#include <stdlib.h>
#include <String.h>
#include <fstream>

#include <iomanip>

#include "DiagInput.h"

#include "TimeTicks.h"

#include "Configured.h"
#include "Compiled.h"
#include "Banner.h"

#include "Verbosity.h"

#include "InternalConfsDiag.h"

#include "FortranLinkage.h"
#include "DiskBuffer.h"

#include "config.h"
#include "VersionDate.h"

using namespace std;

void	MakeBanner()
{
const INT w = 80;
	
	MakeBannerTop(w, "D i a g o n a l i s a t o r");
	center("parallel version based on shared memory", w);
	center(VERSION, w);
	center(DATE, w);
	MakeBannerBottom(w);
}



static MRCIMatrix<double, double>	*mrciMatrix;

extern "C" {


void	FORTRAN_LINKAGE(fmain)(INT *refCSFs, INT *CSFs);

void	FORTRAN_LINKAGE(hmult)(INT *n, double *x, double *y)
{
	cout << "OK, n=" << *n << endl;
INT	SAFs = mrciMatrix->getNExternalsDiag()->getNumberOfTotalSpinAdaptedFunctions();
	cout << "SAFs=" << SAFs << endl;
DiskBuffer	xbuf(SAFs*sizeof(double), "x.dat", 7);
DiskBuffer	ybuf(SAFs*sizeof(double), "y.dat", 7);


	for ( INT i=0 ; i<*n ; i++ )
		xbuf.put(i, x+i*SAFs);

	mrciMatrix->mult(&xbuf, &ybuf, 0, *n-1);

	for ( INT i=0 ; i<*n ; i++ )
		ybuf.get(i, y+i*SAFs);
}


void	FORTRAN_LINKAGE(refmat)(double *h)
{
	mrciMatrix->getRefMatrix(h);
}

void	FORTRAN_LINKAGE(refindex)(INT *i, INT *ind)
{
	*ind = mrciMatrix->getRefIndex(*i-1)+1;
}


};


TimeTicks	globalTime;


void	calc(
	DiagInput & diagInput, INT restart, const char * start, INT slaves)
{
RootEnergies	re;
	re.setStorePTEnergy(diagInput.getStorePTEnergy());
	re.setStorePTCoef(diagInput.getStorePTCoef());
	
Fort31File	f31(diagInput.getMOIntegralFile().name,
				diagInput.getMOIntegralFile().format);
MOIrReps	moirreps(f31);

	cout << "reading configurations..." << endl;
	
ifstream	ConfIn("ConfTree.dat");

	moMapping = MOMapping(ConfIn);

NExternalsDiag	MRDCIConfs(ConfIn, &moirreps);


INT	SAFs = MRDCIConfs.getNumberOfTotalSpinAdaptedFunctions();
INT	RefSAFs = MRDCIConfs.getNumberOfRefConfSpinAdaptedFunctions();


	MRDCIConfs.calcInternalIntersection();


	cout.precision(10);

	{
	ConfigurationSAFNr<MOType> mainA = MRDCIConfs.getConfigurationSAFNr(0);
		if ( diagInput.getStorePTEnergy() && !mainA.getEnergy().getEnergyP() )
		{
			cout << "no PT energies in configuration tree." << endl;
			cout << "be sure you have specified "
				<< "\"StorePTEnergy = yes\" in the selector input." << endl;
			exit(1);
		}
		if ( diagInput.getStorePTCoef() && !mainA.getEnergy().getCICoefP() )
		{
			cout << "no PT coefficients in configuration tree." << endl;
			cout << "be sure you have specified "
				<< "\"StorePTCoef = yes\" in the selector input." << endl;
			exit(1);
		}
	}

	mrciMatrix = new MRCIMatrix<double, double>(&MRDCIConfs,
				f31,
				1024, 1000000,
				slaves);


	FORTRAN_LINKAGE(fmain)(&RefSAFs, &SAFs);


	delete mrciMatrix;
	
}


extern INT writeHamiltonOnly;






int	main(INT argc, char **argv)
{
	globalTime.start();

	MakeBanner();

	StartUp();



INT	slaves = 0;	
INT	readStdin = 0;
INT	restart = 0;
const char	*start = NULL;
	for ( INT i=1 ; i<argc ; i++ )
	{
		if ( argv[i][0]=='-' )
		{
			argv[i][0] = 0;
			switch ( argv[i][1] ) {
			case 'p':
				if ( i+1>=argc )
				{
					cerr << "missing argument to option -p" << endl;
					exit(1);
				}
				slaves = atoi(argv[++i]);
				break;
				
			case 's':
				if ( i+1>=argc )
				{
					cerr << "missing argument to option -s" << endl;
					exit(1);
				}
				start = argv[++i];
				break;

			case 'i':
				readStdin = 1;
				break;
				
			case 'r':
				restart = 1;
				break;

			case 'w':
				writeHamiltonOnly = 1;
				break;

			default:
				cerr << "unknown option " << argv[i] << endl;
				exit(1);
			}
		}
	}
	
/*	if ( slaves )
	{
		cerr << "sorry: this version is non-parallel" << endl;
		exit(1);
	}
*/
	
	if ( slaves )
		cout << "calculation running parallelized on " << slaves << " slaves." << endl;
	else
		cout << "single process mode." << endl;


DiagInput	diagInput(!readStdin);

	if ( readStdin )
	{
		cout << "reading from stdin..." << endl;
		cin >> diagInput;
	}
	
	cout << diagInput;
	cout << endl << endl;
	
	calc(diagInput, restart, start, slaves);

	
}
