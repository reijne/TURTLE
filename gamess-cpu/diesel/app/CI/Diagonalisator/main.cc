#include "../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../lib/QM/MRCIMatrix/MRCIMatrix.h"

#include "../../../lib/QM/MRCIMatrix/DavidsonMatrixStorage.h"

#include "../../../lib/QM/Davidson/DavidsonCI.h"

#include "../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/Configuration/ConfigurationSet.h"

#include "../../common/StartUp.h"
#include "../../../lib/QM/Cache/VarSizeReadOnlyCache.h"


#include <stdlib.h>
#include "../../../lib/Container/String.h"
#include <fstream>

#include <iomanip>

#include "DiagInput.h"

#include "../../../lib/QM/IO/TimeTicks.h"

#include "../../../Configured.h"
#include "Compiled.h"
#include "../../common/Banner.h"
#include "../../..//VersionDate.h"

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../lib/QM/MRTree/Diag/InternalConfsDiag.h"

#include "../../../config.h"

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








TimeTicks	globalTime;
extern INT writeHamiltonOnly;


template <class MatrixType, class VectorType>
void	calc(MRCIMatrix<MatrixType, VectorType>	*mrciMatrix,
	DiagInput & diagInput, INT restart, const char * start, INT slaves)
{
RootEnergies	re;
	re.setStorePTEnergy(diagInput.getStorePTEnergy());
	re.setStorePTCoef(diagInput.getStorePTCoef());

Fort31File	f31(diagInput.getMOIntegralFile());

MOIrReps	moirreps(f31);

	cout << "reading configurations..." << endl;
	
ifstream	ConfIn(diagInput.getConfTreeFileName());
        if (!ConfIn.good()) 
	{
            cout << "Problem opening " << diagInput.getConfTreeFileName() << endl;
            exit(1);
	}

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

	mrciMatrix = new MRCIMatrix<MatrixType, VectorType>(&MRDCIConfs,
				f31,
				1024, 1000000,
				slaves);

//INT	confs = MRDCIConfs.getNumberOfLeaves();


/*	cout << "start." << endl;
	for ( MRTreeIterator i = MRDCIConfs.firstInTree() ;
			!MRDCIConfs.isLastInTree(i) ; MRDCIConfs.nextInTree(i))
		;
	cout << "end." << endl;
*/		
/*	for ( MRTreeIterator i = MRDCIConfs.firstInTree() ;
			!MRDCIConfs.isLastInTree(i) ; MRDCIConfs.nextInTree(i))
	{
//		cout << i.Citer[0].i<< " " << i.Citer[1].i<< " " << i.Citer[2].i<< " " << i.Citer[3].i << ": ";
		cout << MRDCIConfs.getConfigurationSAFNr(i) << endl;
	}
*/	
/*	for ( INT j=0 ; j<MRDCIConfs.getNumberOfLeaves() ; j++ )
	{
		if ( MRDCIConfs.getConfigurationSAFNr(j).isReference() )
		{
			cout << j << ": ";
			cout << MRDCIConfs.getConfigurationSAFNr(j) << endl;
		}
	}
*/
/*	exit(0);
*/
/*	for ( INT i=0 ; i<confs ; i++ )
	{
		ConfigurationSAFNr<MOType> mainA = MRDCIConfs.getConfigurationSAFNr(i);
		cout << mainA << endl;
//		cout << i << ", " << mainA.getSAFNr() << ", " << mainA.getSAFInc() << ":" << mainA << endl;
	}
	exit(0);
*/
/*
Davidson::Type	*p = new Davidson::Type[SAFs*SAFs];

	mrciMatrix->getCIMatrix(p);
Davidson::Type	*pp = p;
	for ( INT i=0 ; i<SAFs ; i++ )
	{
		for ( INT j=0 ; j<SAFs ; j++ )
			cout << *pp++ << " ";
//			printf("%8.4lf ", *pp++);
		cout << endl;
	}

	delete p;
//	exit(0);
*/
/*	{
		CIVectors<Davidson::Type>	*x = new CIVectors<Davidson::Type>(SAFs, 1);
		CIVectors<Davidson::Type>	*y = new CIVectors<Davidson::Type>(SAFs, 1);
		for ( INT i=0 ; i<SAFs ; i++ )
			scanf("%lf", &(*x)(i));

		for ( INT i=0 ; i<SAFs ; i++ )
			printf("%16.10lf\n", (*x)(i));
			
		printf("=================================\n");
		y->clear();
		mrciMatrix->calc(*x, *y);

		for ( INT i=0 ; i<SAFs ; i++ )
			printf("%16.10lf\n", (*y)(i));

		exit(0);
	}
*/




Roots	*roots;

	if ( diagInput.getNumberOfRoots() )
		roots = new Roots(diagInput.getNumberOfRoots(), diagInput.getRootNumbersP());
	else
		roots = new Roots(MRDCIConfs.getNumberOfRoots(), MRDCIConfs.getRootNumbersP());

	cout << endl << endl;
	cout << *roots;
	cout << endl << endl;



DavidsonMatrixStorage<MatrixType, VectorType>	*davidsonMatrixStorage = NULL;
DavidsonCI<MatrixType, VectorType>	*davidson = NULL;




	if ( writeHamiltonOnly )
	{
		davidsonMatrixStorage = new DavidsonMatrixStorage<MatrixType, VectorType>(
			mrciMatrix);

		
		// free integrals
		cout << "freeing integrals." << endl;
		mrciMatrix->freeIntegrals();
		
		delete davidsonMatrixStorage;
		delete mrciMatrix;
		exit(0);
	}




LONG_LONG_INT nonZero = mrciMatrix->estimateNonZeroElements(diagInput.getMaxStorageMem()/
		DavidsonMatrixStorage<MatrixType, VectorType>::getMemory(1)+100);
double	sparsity = (1.0 - 1.0*nonZero/(SAFs*(SAFs+1)/2));
	if ( sparsity>=1.0 )
		cout << "sparsity check aborted" << endl;
	else
		cout << "estimated sparsity: " << 100*sparsity << "%" << endl;

LONG_LONG_INT bytes = DavidsonMatrixStorage<MatrixType, VectorType>::getMemory(nonZero);
	if ( bytes>diagInput.getMaxStorageMem() )
	{
		cout << "available memory (" << (diagInput.getMaxStorageMem() >> 20) <<
			" MB) <= needed memory (" << (bytes >> 20) << " MB) ==>" << endl;
		cout << "using direct method." << endl;

		davidson = new DavidsonCI<MatrixType, VectorType>(SAFs, RefSAFs, *roots,
			mrciMatrix->getNExternalsDiag(), *mrciMatrix, restart,
			diagInput.getRootHoming(),
			((typename DavidsonCI<MatrixType, VectorType>::IterationMode) diagInput.getIterationMode()));
	}
	else
	{
		cout << "available memory (" << (diagInput.getMaxStorageMem() >> 20) <<
			" MB) >= needed memory (" << (bytes >> 20) << " MB) ==>" << endl;
		cout << "storing Hamilton Matrix." << endl;


		davidsonMatrixStorage = new DavidsonMatrixStorage<MatrixType, VectorType>(
			mrciMatrix);

		
		// free integrals
		cout << "freeing integrals." << endl;
		mrciMatrix->freeIntegrals();
		

		davidson = new DavidsonCI<MatrixType, VectorType>(SAFs, RefSAFs, *roots, 
			mrciMatrix->getNExternalsDiag(), 
			*davidsonMatrixStorage, restart,
			diagInput.getRootHoming(),
			((typename DavidsonCI<MatrixType, VectorType>::IterationMode) 
				diagInput.getIterationMode()));

	}


/*
	cout << "checking..." << endl;
MatrixType	*p1 = new MatrixType[SAFs*SAFs];
MatrixType	*p2 = new MatrixType[SAFs*SAFs];
	mrciMatrix->_getTotalMatrix(p1);
	mrciMatrix->getTotalMatrix(p2);
	for ( INT i=0 ; i<SAFs ; i++ )
	{
		for ( INT j=0 ; j<SAFs ; j++ )
			cout << setw(14) << p1[i*SAFs + j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	for ( INT i=0 ; i<SAFs ; i++ )
	{
		for ( INT j=0 ; j<SAFs ; j++ )
			cout << setw(14) << p2[i*SAFs + j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	for ( INT i=0 ; i<SAFs ; i++ )
		for ( INT j=0 ; j<SAFs ; j++ )
		{
			if ( p1[i*SAFs + j] != p2[i*SAFs + j] )
				cout << i << " " << p1[i*SAFs + j]<< " " << p2[i*SAFs + j] << " " << endl;
		}
	cout << "OK" << endl;
*/	





	if ( restart )
	{
		cout << "restarting davidson algorithm" << endl;
		davidson->restart();
	}
	else
		if ( start )
		{
		String	conf(diagInput.getConfTreeFileName());
		String	eigen("Eigenvectors.dat.");
		
			conf += start;
			eigen += start;
			cout << "reading configuration tree and eigenvectors from" << endl;
			cout << conf << " and " << eigen << " respectively."<< endl;
			davidson->start(conf.chars(), eigen.chars());
		}
		else
			davidson->start();


	
	if ( diagInput.getConvergenceEigenvectorChange()!=1 )
		davidson->iterate(
			diagInput.getMaxIters(),
			Davidson<MatrixType, VectorType>::CorrectionVector,
			diagInput.getConvergenceEigenvectorChange());
	else
		davidson->iterate(
			diagInput.getMaxIters(),
			Davidson<MatrixType, VectorType>::Energy,
			diagInput.getConvergenceEnergyChange());
	

VectorType	*RefMatEigenvector[roots->getNumberOfRoots()];
//	cout << roots->getNumberOfRoots() << endl;
	for ( INT ii=0 ; ii<roots->getNumberOfRoots() ; ii++ )
	{
		RefMatEigenvector[ii] = new VectorType[mrciMatrix->getRefDim()];
		memcpy(RefMatEigenvector[ii], 
			davidson->getRefMatEigenvectorP(ii),
			mrciMatrix->getRefDim()*sizeof(VectorType));
	}
		

	delete davidson;
	if ( davidsonMatrixStorage )
		delete davidsonMatrixStorage;



	if ( verbosity.isActive(Verbosity::CacheStatistics) )
		cout << endl << endl << *mrciMatrix->getCache() << endl << endl;
	
	
	delete mrciMatrix;
	
	cout << "=================================================================" << endl;
	cout << endl;
	cout << "                      R e s u l t s" << endl;	
	cout << endl;
	cout << "=================================================================" << endl;
	cout << endl;
	cout << endl;
	cout << *roots;

	cout << endl << endl;
	
	if ( verbosity.isActive(Verbosity::WaveFunction) )
		cout << "wave function with ci^2/nSAFS>" << diagInput.getRefThreshold() << " | refs:" << endl;


Vector<VectorType>	*ev = new Vector<VectorType>(SAFs);
DiskBuffer	*evBuf = new DiskBuffer(SAFs*sizeof(VectorType),
	"Eigenvectors.dat", DiskBuffer::noTempDir);

ConfigurationSet	newRefs;
ConfigurationSet	newPTRefs;

VectorType	EigenvectorOverLap[roots->getNumberOfRoots()*
		roots->getNumberOfRoots()];

	memset(EigenvectorOverLap, 0, roots->getNumberOfRoots()*
		roots->getNumberOfRoots()*sizeof(VectorType));

	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
	{
	EnergyType	RefSum = 0;
	

		if ( verbosity.isActive(Verbosity::WaveFunction) )
		{
			cout << "root # " << i+1 << endl;
			cout << "                ci^2           CSF coef.   status   ";
			for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
				cout << "  PT dE root #" << j+1 << "  ";
			cout << "configuration" << endl;
			cout << "----------------------------------------------------";
			for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
				cout << "------------------";
			cout << "-----------------------------------------" << endl;
		}


		cout.setf(ios::fixed);
		cout.precision(10);

		evBuf->get(i, ev->getP());
INT	safNr = 0;
INT	EigenvectorSAF = 0;
		for ( MRTreeIterator j = MRDCIConfs.firstInTree() ;
				!MRDCIConfs.isLastInTree(j) ; MRDCIConfs.nextInTree(j) )
		{
			ConfigurationSAFNr<MOType> main = MRDCIConfs.getConfigurationSAFNr(j);
		double	c2 = 0;
		double	h;
			for ( INT k=0 ; k<main.getSAFInc() ; k++ )
			{
				h = (*ev)[safNr+k];
				c2 += h*h;
			}

			if ( main.isReference() )
			{
				RefSum += c2;

				for ( INT ii=0 ; ii<roots->getNumberOfRoots() ; ii++ )
					for ( INT k=0 ; k<main.getSAFInc() ; k++ )
					{
/*						cout << EigenvectorOverLap[roots->getNumberOfRoots()*i + ii] << " " <<
							RefMatEigenvector[ii][EigenvectorSAF+k] << " " << (*ev)[safNr+k] << endl;
*/						EigenvectorOverLap[roots->getNumberOfRoots()*i + ii] += 
							RefMatEigenvector[ii][EigenvectorSAF+k]*(*ev)[safNr+k];

					}
				EigenvectorSAF += main.getSAFInc();
			}
			
			if ( verbosity.isActive(Verbosity::WaveFunction) )
			{
//				if ( c2>diagInput.getRefThreshold() || main.isReference() )
//				changed: 16.01.1998
				if ( c2/main.getSAFInc()>diagInput.getRefThreshold() || main.isReference() )
				{
					cout << setw(20) << c2 <<
						setw(20) << (*ev)[safNr+0] << ": " 
							<< main;
					if ( main.isReference() )
					{
					InternalConfsDiag *internal = MRDCIConfs.operator [] (0);
					ContainerIterator iteri = internal->first();
					INT	n = 1;
						while ( (!(*internal)[iteri]->isReference()) || !(main == *(*internal)[iteri]) )
						{
							if ( (*internal)[iteri]->isReference() )
								n++;
							internal->next(iteri);
						}
						cout << "\t\t# " << n;
					}
					cout << endl;
					for ( INT k=1 ; k<main.getSAFInc() ; k++ )
					cout << setw(20) << " " <<
						setw(20) << (*ev)[safNr+k] << endl;
				}
			}
			
//			if ( c2>diagInput.getRefThreshold() )
//			changed: 16.01.1998
			if ( c2/main.getSAFInc()>diagInput.getRefThreshold() )
					newRefs.add(main);
				
//			if ( c2>diagInput.getPTRefThreshold() )
//			changed: 16.01.1998
			if ( c2/main.getSAFInc()>diagInput.getPTRefThreshold() )
					newPTRefs.add(main);
				
			safNr += main.getSAFInc();
		}
		if ( verbosity.isActive(Verbosity::WaveFunction) )
		{
			cout << endl;
			cout << "ci^2 of references: " << RefSum << endl;
			cout << endl;
			cout << "===========================================" << endl;
			cout << endl;
		}
	}
	delete ev;
	delete evBuf;

	for ( INT ii=0 ; ii<roots->getNumberOfRoots() ; ii++ )
		delete[] RefMatEigenvector[ii];
		

	cout << "eigenvector overlap with reference space:" << endl;
	cout << "(column: reference root, row: MR-CI root)" << endl;
	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
	{
		for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
			cout << setw(14) << EigenvectorOverLap[roots->getNumberOfRoots()*i + j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << endl;
INT	flipped = 0;
	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
	{
	VectorType	sum = 0;
		
		for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
			sum += EigenvectorOverLap[roots->getNumberOfRoots()*i + j]*
				EigenvectorOverLap[roots->getNumberOfRoots()*i + j];

		if ( fabs(EigenvectorOverLap[roots->getNumberOfRoots()*i + i])<0.8 )
			flipped = 1;
		if ( fabs(sum)<0.1 )
			cout << "***** warning: small eigenvector overlap for root #" << i << " *****" << endl;
		else
		{
			sum = sqrt(1/sum);
			for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
				EigenvectorOverLap[roots->getNumberOfRoots()*i + j] *= sum;
		}
	}
	cout << "renormalized eigenvector overlap with reference space:" << endl;
	cout << "(column: reference root, row: MR-CI root)" << endl;
	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
	{
		for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
			cout << setw(14) << EigenvectorOverLap[roots->getNumberOfRoots()*i + j] << " ";
		cout << endl;
	}
	cout << endl;
	cout << endl;
	if ( flipped )
	{
		cout << "********************************************************************" << endl;
		cout << endl;
		cout << "                attention: probably flipped roots" << endl;
		cout << endl;
		cout << "********************************************************************" << endl;
		cout << endl;
		cout << endl;
	}
	
	cout << "configurations above reference threshold (" << diagInput.getRefThreshold() << "):" << endl;
	cout << "{" << endl;
	newRefs.writeToStream(cout);
	cout << "}" << endl;
	cout << endl;
	cout << endl;
	cout << "configurations above PT reference threshold (" << diagInput.getPTRefThreshold() << "):" << endl;
	cout << "{" << endl;
	newPTRefs.writeToStream(cout);
	cout << "}" << endl;
}




int	main(int argc, char **argv)
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
	
	switch ( diagInput.getPrecision() ) {
	case DiagInput::floatPrec:
		{
		typedef float VectorType;
		typedef float MatrixType;

		MRCIMatrix<MatrixType, VectorType>	*mrciMatrix = NULL;
			calc(mrciMatrix, diagInput, restart, start, slaves);
		}
		break;
		
	case DiagInput::doublePrec:
		{
		typedef double VectorType;
		typedef double MatrixType;

		MRCIMatrix<MatrixType, VectorType>	*mrciMatrix = NULL;
			calc(mrciMatrix, diagInput, restart, start, slaves);
		}
		break;

	case DiagInput::undefinedPrec:
		break;
	}
	
}
