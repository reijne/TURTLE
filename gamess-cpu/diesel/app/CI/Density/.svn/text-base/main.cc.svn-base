#include "../../../lib/QM/MRTree/Diag/NExternalsDiag.h"


#include "../../../lib/QM/Davidson/Davidson.h"

#include "../../../lib/QM/MRTree/MRTreeIterator.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/MRCIMatrix/MRCIMatrix.h"
#include "../../../lib/QM/RepresentationMatrices/CIVectors.h"

#include "../../../lib/QM/IntegralContainer/DensityMatrix.h"

#include "../../../lib/Container/GrepAwk.h"


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



#include <stdlib.h>
#include "../../../lib/Container/String.h"
#include <fstream>

#include <iomanip>


#include "../../../config.h"
#include "../../../VersionDate.h"

using namespace std;

int	main(int argc, char **argv)
{
	cout << "dens (Part of DIESEL-MR-CI), " << VERSION << ", " << DATE << endl << endl;
	
	if ( argc<6 || argc>8 )
	{
		cerr << endl;
		cerr << "purpose:" << endl;
		cerr << "the program calculates the one electron density matrix" << endl; 
		cerr << endl;
		cerr << "usage: dens [-n -o -t -r] motra-input-file Thresh {#Lstate1[-#LstateN], all} Ldir {[#Rstate1[-RstateM], all}  [Rdir]]" << endl;
		cerr << endl;
		cerr << endl;
		exit(1);
	}

	
INT	stateLs = 0;
INT	stateLe = 0;
INT	stateRs = 0;
INT	stateRe = 0;
int 	stateLs_2,stateLe_2;
int	stateRs_2,stateRe_2;

	if ( !strcmp(argv[4], "all") )
		stateLe = stateLs = -1;
	else
		if ( sscanf(argv[4], "%d-%d", &stateLs_2, &stateLe_2)!=2 ) {
			stateLe = stateLs_2;
		        stateLs=stateLs_2;}
 		else {
			stateLe=stateLe_2;
			stateLs=stateLs_2;
		}
		

	if ( argc>=6 )
	{
		if ( !strcmp(argv[6], "all") )
			stateRe = stateRs = -1;
		else
			if ( sscanf(argv[6], "%d-%d", &stateRs_2, &stateRe_2)!=2 ){
				stateRs = stateRs_2;
				stateRe = stateRs_2;
			} else {
				stateRs = stateRs_2;
				stateRe = stateRe_2;
			}
	}

	{
	ifstream	f31("fort.31");
		if ( !f31 )
		{
			cerr << "no file \"fort.31\"" << endl;
			exit(1);
		}
		if ( argv[1][0]!='-' )
		{
			cerr << "file format argument missing." << endl;
			exit(1);
		}
	}
	
String	motraInputFilename(argv[2]);

String	Thresh(argv[3]);
	
String	LdirFilename(argv[5]);
String	LConfTreeFilename(LdirFilename+"/ConfTree.dat."+Thresh);
String	LEigenvectorsFilename(LdirFilename+"/Eigenvectors.dat."+Thresh);



String	RdirFilename(LdirFilename);
String	RConfTreeFilename(LConfTreeFilename);
String	REigenvectorsFilename(LEigenvectorsFilename);

INT	readSecondConfTree = argc==8;
INT	readSecondEigenvectors = argc>=7;

	if ( readSecondConfTree )
	{	
		RdirFilename = argv[7];
		RConfTreeFilename = RdirFilename+"/ConfTree.dat."+Thresh;
		REigenvectorsFilename = RdirFilename+"/Eigenvectors.dat."+Thresh;
	}
	
	
	
Fort31File	*f31 = NULL;
	switch ( argv[1][1] )
	{
	case 'n':
		f31 = new Fort31File("fort.31", Fort31RecordFormatNew);
		break;

	case 'o':
		f31 = new Fort31File("fort.31", Fort31RecordFormatOld);
		break;

	case 't':
		f31 = new Fort31File("fort.31", Fort31RecordFormatTRADPT);
		break;
	
	case 'r':
		f31 = new Fort31File("fort.31", RIFormat);
		break;
	
	default:
		cerr << "wrong file format" << endl;
		exit(1);
	}

MOIrReps	moirreps(*f31);

INT* nFrozen = new INT[moirreps.getNumberOfIrReps()];
INT* nDeleted = new INT[moirreps.getNumberOfIrReps()];


if (f31->format != RIFormat)
{
ifstream	motraInput(motraInputFilename);	
		if ( !motraInput )
		{
			cerr << "no file \"" << motraInputFilename << "\"" << endl;
			exit(1);
		}

GrepAwk	ga(motraInput);

	ga.setIgnoreCase(1);

	if ( !ga.grep("&MOTRA") )
	{
		cerr << "motra.in: &MOTRA not found\n" << endl;
		exit(1);
	}
	if ( !ga.grep("frozen") )
	{
		cerr << "motra.in: frozen not found\n" << endl;
		exit(1);
	}
	ga++;
	for ( INT i=0 ; i<moirreps.getNumberOfIrReps() ; i++ )
		nFrozen[i] = atoi(ga.getWord(i+1));

	if ( !ga.grep("delete") )
	{
		cerr << "motra.in: delete not found\n" << endl;
		exit(1);
	}
	ga++;
	for ( INT i=0 ; i<moirreps.getNumberOfIrReps() ; i++ )
		nDeleted[i] = atoi(ga.getWord(i+1));
	
}
else
{
	//cout << "Die Infos MOs -frozen/deleted fehlen " << endl;
	//cout << "fuer RI-Rechnung ist das aber vollkommen ausreichend, weil Infos in files coredata und cisdata" << endl << endl;
	cout << "please recognize, that ridens does not write any infos" << endl;
	cout << "about frozen and deleted orbitals." << endl;
	cout << "These Informations are included in some other files and" << endl;
	cout << "are taken upon request." << endl;
	// In einer sp"ateren Version w"are eine Erweiterung denkbar. Die n"otigen Infos k"onnen mit Hilfe der
	// im Header risym.h befindlichen Funktionen eingelesen werden.
	for (INT i=0; i<moirreps.getNumberOfIrReps() ; i++ )
	{
		nFrozen[i] = 0;
		nDeleted[i] = 0;
	}
}

	cout << endl;
	cout << "reading configurations for first state..." << endl;
	
ifstream	ConfInL(LConfTreeFilename);
MOMapping moMappingDummy(ConfInL);

NExternalsDiag	*MRDCIConfsL = new NExternalsDiag(ConfInL, &moirreps);


INT	SAFsL = MRDCIConfsL->getNumberOfTotalSpinAdaptedFunctions();
INT	SAFsR = SAFsL;


typedef	double	VectorType;
typedef double	MatrixType;

DiskBuffer	*evBufL = new DiskBuffer(SAFsL*sizeof(VectorType),
			LEigenvectorsFilename, DiskBuffer::noTempDir);
	if ( stateLs==-1 )
	{
		stateLs = 1;
		stateLe = evBufL->getNumberOfObjects();
	}
	cout << ":::" << SAFsL << " " << evBufL->getNumberOfObjects() << endl;
CIVectors<VectorType>	*evL = new CIVectors<VectorType>(SAFsL, stateLe-stateLs+1);

	for ( INT state=stateLs, istate=0 ; state<=stateLe ; state++, istate++ )
	{
		cout << "reading eigenvector for first state (" << state << ")...";
		evBufL->get(state-1, evL->getP(istate));

	}

MRMOs	mrmosL = *MRDCIConfsL->getMRMOs();
	moMappingDummy.setInternal(&mrmosL);
	mrmosL.initIntExt();



MRMOs	mrmosR;

NExternalsDiag	*MRDCIConfsR = MRDCIConfsL;
CIVectors<VectorType>	*evR = NULL;
DiskBuffer	*evBufR = NULL;

	cout << endl << endl;

	if ( readSecondConfTree )
	{
		cout << endl;
		cout << "reading configurations for second state..." << endl;

	ifstream	ConfInR(RConfTreeFilename);
		MOMapping moMappingDummy(ConfInR);
		MRDCIConfsR = new NExternalsDiag(ConfInR, &moirreps);

		SAFsR = MRDCIConfsR->getNumberOfTotalSpinAdaptedFunctions();

		mrmosR = *MRDCIConfsR->getMRMOs();
		moMappingDummy.setInternal(&mrmosR);

		cout << endl << endl;
	}

	if ( readSecondEigenvectors )
	{
		evBufR = new DiskBuffer(SAFsR*sizeof(VectorType), 
			REigenvectorsFilename, DiskBuffer::noTempDir);
		if ( stateRs==-1 )
		{
			stateRs = 1;
			stateRe = evBufR->getNumberOfObjects();
		}
		evR = new CIVectors<VectorType>(SAFsR, stateRe-stateRs+1);

		for ( INT state=stateRs, istate=0 ; state<=stateRe ; state++ , istate++ )
		{
			cout << "reading eigenvector for second state (" << state << ")...";
			evBufR->get(state-1, evR->getP(istate));

		}

		cout << endl << endl;
	}

	if ( readSecondConfTree )
	{
		for ( INT i=0 ; i<mrmosL.getMaxMO() ; i++ )
			if ( mrmosR.isInternal(i+1) )
				mrmosL.setInternal(i+1);

		mrmosL.initIntExt();

		for ( INT i=0 ; i<mrmosL.getMaxMO() ; i++ )
		{
			if ( mrmosL.isInternal(i+1)!=MRDCIConfsR->getMRMOs()->isInternal(i+1) )
			{
				cout << "reorganizing internal space of second wave function..." << flush;
				moMapping.reorderNoHoles(MRDCIConfsR, evR, &mrmosL);
				break;
			}
		}


		for ( INT i=0 ; i<mrmosL.getMaxMO() ; i++ )
		{
//			cout << i+1 << " " << mrmosL.isInternal(i+1) << " " << MRDCIConfsL->getMRMOs()->isInternal(i+1) << endl;
			if ( mrmosL.isInternal(i+1)!=MRDCIConfsL->getMRMOs()->isInternal(i+1) )
			{
				cout << "reorganizing internal space of first wave function..." << flush;
				moMapping.reorderNoHoles(MRDCIConfsL, evL, &mrmosL);
				break;
			}
		}
	}



	cout.precision(10);
	
	

INT maxMO = moirreps.getMaxMO();
DensityMatrix::Type	**densMat = NULL;

	if ( evR )
		densMat = new DensityMatrix::Type * [evL->getN()*evR->getN()];
	else
		densMat = new DensityMatrix::Type * [evL->getN()];
		
	{
	INT	nDens = 0;
		for ( INT ix=0 ; ix<evL->getN() ; ix++ )
		{
		INT	iys = ix;
		INT	iye = ix+1;
			if ( evR )
			{
				iys = 0;
				iye = evR->getN();
			}
			for ( INT iy=iys ; iy<iye ; iy++ )
				densMat[nDens++] = new DensityMatrix::Type[maxMO*maxMO];
		}
	}


MRCIMatrix<MatrixType, VectorType>	*mrciMatrix = new 
	MRCIMatrix<MatrixType, VectorType>(MRDCIConfsL, MRDCIConfsR, *f31);

	cout << endl;
	cout << "calculating density matrix(ces):" << endl;

	mrciMatrix->dens(evL, evR, densMat);

/*	for ( INT ix=0 ; ix<evL->getN() ; ix++ )
	{
		for ( INT i=0 ; i<maxMO ; i++ )
		{
			for ( INT j=0 ; j<maxMO ; j++ )
			{
				cout << setprecision(8) << setw(16) << densMat[ix][i*maxMO+j];
			}
			cout << endl;
		}

		cout << "---------------" << endl;
	}
*/
	mkdir("Densities", S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
	{
	INT	nDens = 0;
		for ( INT ix=0 ; ix<evL->getN() ; ix++ )
		{
		INT	iys = ix;
		INT	iye = ix+1;
			if ( readSecondEigenvectors  )
			{
				iys = 0;
				iye = evR->getN();
			}
			if ( (!readSecondConfTree) && readSecondEigenvectors )
				nDens += ix;
			for ( INT iy=((!readSecondConfTree) && readSecondEigenvectors ? ix : iys) ; iy<iye ; iy++ )
			{
			DensityMatrix::Type	*densMatReal = new DensityMatrix::Type[maxMO*maxMO];


				for ( INT i=0 ; i<maxMO ; i++ )
					for ( INT j=0 ; j<maxMO ; j++ )
						densMatReal[(moMapping.getReal(i+1)-1)*maxMO+moMapping.getReal(j+1)-1] = 
							densMat[nDens][i*maxMO+j];


				delete densMat[nDens++];


			DensityMatrix	densityMatrix(moirreps, nFrozen, nDeleted, densMatReal,
					readSecondConfTree || readSecondEigenvectors && ix!=iy);


			double	sum = 0;
				for ( INT i=0 ; i<maxMO ; i++ )
						sum += densMatReal[i*maxMO+i];
				cout << "trace=" << sum << endl;


				delete densMatReal;


			char	fname[100];
				sprintf(fname, "Densities/Density.dat.I%sR%d_I%sR%d.%s", 
					LdirFilename.chars(), ix+stateLs, RdirFilename.chars(), iy+stateLs,
					Thresh.chars());

//	cout << densityMatrix << endl;
			ofstream	f(fname);
				densityMatrix.writeToStream(f);
			}

		}
	}
	delete densMat;
	delete[] nFrozen;
	delete[] nDeleted;


//	cout << densityMatrix << endl;



/*	for ( INT i=0 ; i<maxMO ; i++ )
	{
		for ( INT j=0 ; j<maxMO ; j++ )
		{
			cout << setprecision(8) << setw(16) << densMatReal[i*maxMO+j];
		}
		cout << endl;
	}

	cout << "---------------" << endl;
	cout << densityMatrix << endl;
*/

//********************************************************************
//********************************************************************
//********************************************************************
//	for checking:
/*
ifstream	s("CIDENS");
const INT maxchars = 200;
char	buf[maxchars];

char	cmpstring[1000];
	sprintf(cmpstring, "CI DENSITY MATRIX FOR STATE NO.     2");
	while ( s.getline(buf, maxchars) )
	{
		if ( !strcmp(buf, cmpstring) )
			break;
	}
	
	s.getline(buf, maxchars);
	s.getline(buf, maxchars);
	s.getline(buf, maxchars);
	s.getline(buf, maxchars);

double	*densMatOldCI = new double[maxMO*maxMO];
double	*p = densMatOldCI;

	cout << "MaxMO=" << maxMO << endl;

	for ( INT k=0 ; k<maxMO*maxMO ; k++ )
	{
		s.read(buf, 22);
		buf[22] = 0;
		*strchr(buf, 'D') = 'E';
		*p++ = atof(buf);
		if ( k%4==3 || k==maxMO*maxMO-1 )
			s.getline(buf, maxchars);
	}
	

	for ( INT i=0 ; i<maxMO ; i++ )
	{
		for ( INT j=0 ; j<maxMO ; j++ )
		{
			if ( fabs(fabs(densMatReal[i*maxMO+j])-fabs(densMatOldCI[i*maxMO+j]))>1e-5 )
			cout << "(" << i+1 << "," << j+1 << ")" <<
			 setprecision(8) << setw(16) << densMatReal[i*maxMO+j] <<
			 setprecision(8) << setw(16) << densMatOldCI[i*maxMO+j] << endl;
		}
		cout << endl;
	}
*/
//********************************************************************
//********************************************************************
//********************************************************************



	delete mrciMatrix;
	delete evL;
	delete evBufL;
	delete MRDCIConfsL;
	if ( readSecondEigenvectors )
	{
		delete evR;
		delete evBufR;
	}
	if ( readSecondConfTree )
		delete MRDCIConfsR;
		
	delete f31;
}
