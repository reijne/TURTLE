//***********************************************************************
//
//	Name:			DiagIntContainer.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************



#include "DiagIntContainer.h"

#include "../IO/Fortran/FortranFileIO.h"
#include "../IO/Fortran/Fort31FirstRecord.h"
#include "../IO/Fortran/Fort31SecondRecord.h"



#include "../RepresentationMatrices/HMatElements.h"
#include "../Configuration/TableCase.h"
#include "../Cache/VarSizeReadOnlyCache.h"
#include "../RepresentationMatrices/RepresentationMatrixFortranInterface.h"


#include "../IO/TimeTicks.h"
#include "../MO/MOType.h"
#include "../MO/Iterators/MOEquivalence.h"


#include <iomanip>

#include <stdlib.h>
#include "RIFourIndexIntegralContainer.h"
#include "RITwoIndexIntegralContainer.h"

using namespace std;


DiagIntContainer::DiagIntContainer(Fort31File _f31, INT Multiplicity)
{
	f31 = _f31;
	core = 0;
	moirreps = new MOIrReps(_f31);

	maxMO = moirreps->getMaxMO();
	
	active = new INT[maxMO+1];
	memset(active, 0, (maxMO+1)*sizeof(INT));
	
	nIrreps = moirreps->getNumberOfIrReps();

	pOne = new OneElectronIntegralType[maxMO];
	pTwoIIII = new TwoElectronIntegralType[maxMO];
	pTwoIIJJ = new TwoElectronIntegralType[maxMO*maxMO];
	pTwoIJJI = new TwoElectronIntegralType[maxMO*maxMO];
	if(f31.format != RIFormat)
	{
		load1Integrals();
		load2Integrals();
	}
	else
	{
		//	load1Integrals
		cout << "read RI-Information to construct one-electron integrals" << endl;
		MRMOs	mrmos(f31);
		RITwoIndexIntegralContainer	ri2c(mrmos);
		ri2c.loadIntegrals(f31);
		OneElectronIntegralIndex	e1Index(0,0);
		OneElectronIntegralType		e1IntegralValue= 0;
		INT	intsContained = 0;
		for (INT i=1; i<=maxMO; i++)
		{
			e1Index.set(i,i);
			e1IntegralValue = ri2c.operator[](e1Index);
			setOneInt(i,e1IntegralValue);
			intsContained++;
		}
		cout << "number of read one electron integrals: " << intsContained << endl;
		//	read core
		core = ri2c.getCore();
		
		//	load2Integrals
		cout << "read RI-Information to construct two-electron integrals" << endl;
		
		RIFourIndexIntegralContainer	ri4c(mrmos);
		//RIFourIndexIntegralContainer	ri4c(mrmos,100);
		//#include "SharedMemory.h"
		//INT	size = ri4c.getDataSize();
		//Int4sharedMem = new SharedMemory(size, SharedMemory::StandAlone);
		//ri4c.setSharedMem(Int4sharedMem);
		//ri4c.allocate();
		ri4c.loadIntegrals(f31);
		TwoElectronIntegralIndex<MOType> 	Index(0,0,0,0);
		IntegralType	IntegralValue = 0;
		intsContained = 0;
		for (INT i=1; i<=maxMO; i++)
		{
			Index.set(i,i,i,i);
			ri4c.get(Index,IntegralValue);
			this->setTwoIntIIII(i, IntegralValue);
			intsContained++;
			
			for (INT j=1; j<=i; j++)
			{
				Index.set(i,i,j,j);
				ri4c.get(Index,IntegralValue);			
				this->setTwoIntIIJJ(i,j, IntegralValue);
				intsContained++;

				Index.set(i,j,i,j);
				ri4c.get(Index,IntegralValue);			
				this->setTwoIntIJJI(i,j, IntegralValue);
				intsContained++;
			}
		}
		cout << "number of stored two electron integrals: " << intsContained << endl;
	}
	




////////////////////////////////////////////////////////////////////////
//
//	initialize representation matrices for diagonal case
//

//HMatElements init1(intContainer, int2Container, core);

RepresentationMatrixFortranInterface<MatrixType>	repMatFInt(Multiplicity);

	repMatsP5 = new HMatElements<MatrixType, VectorType> * [MAXSGATABLEOPENSHELLS+1];
	memset(repMatsP5, 0, (MAXSGATABLEOPENSHELLS+1)*sizeof(HMatElements<MatrixType, VectorType> *));
	
BinomialCoefficient	binom(MAXOPENSHELLS);
TableCase<MOType>	tablecase(&binom);				// table case
VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	cache(
			1 << 10, 1 << 20);

	minOpenShells = Multiplicity-1;
	maxOpenShells = MAXOPENSHELLS;

	for ( INT i=minOpenShells ; i<=MAXSGATABLEOPENSHELLS ; i+=2 )
	{
//		cout << i << " open shells:" << endl;
		tablecase.setNr(i, 28);	// P=5, case 14
		repMatsP5[i] =  new HMatElements <MatrixType, VectorType>
			(*(cache)[TableKey(tablecase)]);
//		cout << *repMatsP5[i] << endl;
	}

}


DiagIntContainer::~DiagIntContainer()
{
	delete active;
	delete pOne;
	delete pTwoIIII;
	delete pTwoIIJJ;
	delete pTwoIJJI;
	delete moirreps;

	for ( INT i=0 ; i<MAXSGATABLEOPENSHELLS ; i++ )
		if ( repMatsP5[i] )
			delete repMatsP5[i];
	delete repMatsP5;
}



void	DiagIntContainer::load2Integrals()
{	
////////////////////////////////////////////////////////////////////////////
//
//	two electron integrals	
//	


FortranFileIO	fort31(f31.name);
INT	integrals = 0, intsContained = 0;

TimeTicks	ticks;
	ticks.start();

MOType	ind[4];

	switch ( f31.format ) {
	case GAMESSC1:
		{
			fort31.nextRecord();
			fort31.nextRecord();
			fort31.nextRecord();
		const INT recsize = 15000;
		struct Record {
			INT n; 
			unsigned char label[recsize][4];
			double	a[recsize];
			} data;

		INT	startMO = 0;
			{
				fort31.read(&data, sizeof(Record));
				startMO = data.label[0][0] -1;
				fort31.previousRecord();
			}


			while ( 1 ) 
			{
				fort31.read(&data, sizeof(Record));
				data.n = -data.n;
				if ( data.n<1 )
					break;

				for ( INT i=0 ; i<data.n ; ++i )
				{
					ind[0] = moMapping.getContinuous(data.label[i][3] - startMO);
					ind[1] = moMapping.getContinuous(data.label[i][2] - startMO);
					ind[2] = moMapping.getContinuous(data.label[i][1] - startMO);
					ind[3] = moMapping.getContinuous(data.label[i][0] - startMO);
					intsContained += this->set(ind, data.a[i]);
					
					integrals++;
				}
			}
		}
		break;
	default:
	{

INT	recSize = 2000;
IntegralType* buf = new IntegralType[recSize];


Fort31SecondRecord	f31Second(f31.name, f31.format);
TimeTicks	ticks;

	intsContained = 0;
	cout << "reading two electron integrals...";
	cout.flush();
	ticks.start();
	fort31.rewind();
	if ( f31.format==Fort31RecordFormatTRADPT )
		fort31.nextRecord();
	fort31.nextRecord();
	fort31.nextRecord();
	integrals = 0;

INT	pos = recSize;

	for ( INT i=0 ; i<moirreps->getNumberOfIrReps() ; i++ )
		for ( INT j=0 ; j<=i ; j++ )
			for ( INT k=0 ; k<=i ; k++ )
				for ( INT l=0 ; l<=((i==k) ? j : k) ; l++ )
				{//	cout << "(" << i << " " << j << "|" << k << " " << l << ")" << endl;
										
					
					if ( moirreps->getProd(
							moirreps->getProd(i, j), moirreps->getProd(k, l)) )
						continue;
					if ( i==l )	
					{
//					cout << "(AA|AA)" << endl;
	// (AA|AA) ---------------------------------------------------------
	for ( INT ii=0 ; ii<moirreps->getInIrRep(i) ; ii++ )
	{	ind[0] = ii + moirreps->getStartMO(i);
		for ( INT jj=0 ; jj<=ii ; jj++ )
		{	ind[1] = jj + moirreps->getStartMO(i);
			for ( INT kk=0 ; kk<=ii ; kk++ )
			{	ind[2] = kk + moirreps->getStartMO(i);
				for ( INT ll=0 ; ll<=((ii==kk) ? jj : kk) ; ll++ )
				{	ind[3] = ll + moirreps->getStartMO(i);
					
					if ( pos>=recSize )
					{	pos = 0;
						fort31.read(buf, recSize*sizeof(IntegralType));
					}
					intsContained += this->set(ind, buf[pos++]);
					
					integrals++;
				}
			}
		}
	}
					}
					else
					if ( i==j )
					{
//					cout << "(AA|BB)" << endl;
	// (AA|BB) ---------------------------------------------------------
	for ( INT ii=0 ; ii<moirreps->getInIrRep(i) ; ii++ )
	{	ind[0] = ii + moirreps->getStartMO(i);
		for ( INT jj=0 ; jj<=ii ; jj++ )
		{	ind[1] = jj + moirreps->getStartMO(i);
			for ( INT kk=0 ; kk<moirreps->getInIrRep(k) ; kk++ )
			{	ind[2] = kk + moirreps->getStartMO(k);
				for ( INT ll=0 ; ll<=kk ; ll++ )
				{	ind[3] = ll + moirreps->getStartMO(k);
					
					if ( pos>=recSize )
					{	pos = 0;
						fort31.read(buf, recSize*sizeof(IntegralType));
					}
					
					intsContained += this->set(ind, buf[pos++]);
					integrals++;
				}
			}
		}
	}
					}
					else
					if ( i==k )
					{
//					cout << "(AB|AB)" << endl;
	// (AB|AB) ---------------------------------------------------------
	for ( INT ii=0 ; ii<moirreps->getInIrRep(i) ; ii++ )
	{	ind[0] = ii + moirreps->getStartMO(i);
		for ( INT jj=0 ; jj<moirreps->getInIrRep(j) ; jj++ )
		{	ind[1] = jj + moirreps->getStartMO(j);
			for ( INT kk=0 ; kk<=ii ; kk++ )
			{	ind[2] = kk + moirreps->getStartMO(i);
				for ( INT ll=0 ; ll<moirreps->getInIrRep(j) ; ll++ )
				{	if ( kk==ii && ll>jj )
						continue;

					ind[3] = ll + moirreps->getStartMO(j);

					if ( pos>=recSize )
					{	pos = 0;
						fort31.read(buf, recSize*sizeof(IntegralType));
					}
					
					intsContained += this->set(ind, buf[pos++]);
					integrals++;
				}
			}
		}
	}
					}
					else
					{
//					cout << "(AB|CD)" << endl;
	// (AB|CD) ---------------------------------------------------------
	for ( INT ii=0 ; ii<moirreps->getInIrRep(i) ; ii++ )
	{	ind[0] = ii + moirreps->getStartMO(i);
		for ( INT jj=0 ; jj<moirreps->getInIrRep(j) ; jj++ )
		{	ind[1] = jj + moirreps->getStartMO(j);
			for ( INT kk=0 ; kk<moirreps->getInIrRep(k) ; kk++ )
			{	ind[2] = kk + moirreps->getStartMO(k);
				for ( INT ll=0 ; ll<moirreps->getInIrRep(l) ; ll++ )
				{	ind[3] = ll + moirreps->getStartMO(l);

					if ( pos>=recSize )
					{	pos = 0;
						fort31.read(buf, recSize*sizeof(IntegralType));
					}
					pos++;
					integrals++;
				}
			}
		}
	}

					} // end if
	}	
	delete[] buf;
	}
	}
	ticks.stop();
	cout << ticks << endl;
	cout << "number of read two electron integrals: " << integrals << endl;
	cout << "number of stored two electron integrals: " << intsContained << endl;
}	
	


void	DiagIntContainer::load1Integrals()
{	
////////////////////////////////////////////////////////////////////////////
//
//	one electron integrals	
//	
	
FortranFileIO	fort31(f31.name);
		
INT	integrals = 0;

TimeTicks	ticks;

	cout << "reading one electron integrals...";
	cout.flush();
	ticks.start();
	switch ( f31.format ) {
	case GAMESSC1:
		{
			fort31.nextRecord();
			fort31.read(&core);
			
		INT	recSize = moirreps->getMaxMO()*(moirreps->getMaxMO()+1)/2;
		OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];

		OneElectronIntegralType	*p = buf;

			fort31.read(buf, recSize*sizeof(OneElectronIntegralType));


			for ( INT i=0 ; i<moirreps->getMaxMO() ; i++ )
			{	
				for ( INT j=0 ; j<=i ; j++ )
				{
					if ( moirreps->getProd(moirreps->getIrRep(i+1), moirreps->getIrRep(j+1)) )
						continue;

					if ( i==j )
						setOneInt(moMapping.getContinuous(i+1), *p);

					integrals++;
					p++;
				}
			}
		delete[] buf;
		}
		break;

	case Fort31RecordFormatNew:
		{
		Fort31FirstRecord	f31_1(f31.name, f31.format);
		Fort31SecondRecord	f31Second(f31.name, f31.format);
		INT	recSize = f31_1.getN1el();
		OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];

			fort31.end();
			fort31.previousRecord();
			fort31.previousRecord();
			fort31.previousRecord();
			fort31.previousRecord();
			integrals = 0;

			fort31.read(buf, recSize*sizeof(OneElectronIntegralType));



		OneElectronIntegralType	*p = buf;
			for ( INT i=0 ; i<moirreps->getMaxMO() ; i++ )
			{
				for ( INT j=0 ; j<=i ; j++ )
				{	if ( moirreps->getProd(moirreps->getIrRep(i+1), moirreps->getIrRep(j+1)) )
						continue;

					if ( i==j )
						setOneInt(i+1, *p);

					integrals++;
					p++;
				}
			}

			fort31.end();
			fort31.previousRecord();

			fort31.read(&core, sizeof(double));
			delete[] buf;
		}
		break;

	case Fort31RecordFormatTRADPT:
		{
		Fort31FirstRecord	f31_1(f31.name, f31.format);
		Fort31SecondRecord	f31Second(f31.name, f31.format);
		INT	recSize = 2000;
		OneElectronIntegralType* buf = new OneElectronIntegralType[recSize];
		INT	pos = recSize;

			fort31.end();
			// skip VNUC
			fort31.previousRecord();

		INT	l = f31_1.getNBOX()*(f31_1.getNBOX()+1)/2;
			for ( INT i=0 ; i<(l+recSize-1)/recSize ; i++ )
			{
				fort31.previousRecord();
				fort31.previousRecord();
			}


			for ( INT i=0 ; i<moirreps->getMaxMO() ; i++ )
			{
				for ( INT j=0 ; j<=i ; j++ )
				{	if ( moirreps->getProd(moirreps->getIrRep(i+1), moirreps->getIrRep(j+1)) )
						continue;

					if ( pos>=recSize )
					{	pos = 0;
						fort31.read(buf, recSize*sizeof(OneElectronIntegralType));
					}

					if ( i==j )
						setOneInt(i+1, buf[pos]);
					pos++;
					integrals++;
				}
			}

		fort31.end();
		fort31.previousRecord();
		fort31.read(&core, sizeof(double));
		delete[] buf;
		}
		break;
		
	default:
		exit(1);

	}


	ticks.stop();
	cout << ticks << endl;

	cout << "number of read one electron integrals: " << integrals << endl;
}



EnergyType	DiagIntContainer::calcDiagClosed(INT *occupation)
{
/*
	nMOs = 0;
	for ( INT irrep=0 ; irrep<moirreps->getNumberOfIrReps() ; irrep++ )
		nMOs += occupation[irrep];
*/	
MOType* mos = new MOType[nMOs];

INT	i = 0;
	for ( INT irrep=0 ; irrep<moirreps->getNumberOfIrReps() ; irrep++ )
	{
	// optimization bug in GCC 2.8.0 !!!!!
	MOType	m = occupation[irrep];
		for ( MOType mo=0 ; mo<m ; mo++ )
			mos[i++] = mo + moirreps->getStartMO(irrep);
	}
	return calcDiagClosed(mos, nMOs);
	delete[] mos;
}


EnergyType	DiagIntContainer::calcDiagClosed(const MOType *mos, INT nMOs) const
{
EnergyType	diag = 0;

	for ( INT i=0 ; i<nMOs ; i++ )
	{
	MOType	mo1 = mos[i];
		diag += 2*getOneInt(mo1) + getTwoIntIIII(mo1);
		for ( INT j=i+1 ; j<nMOs ; j++ )
		{
		MOType	mo2 = mos[j];
			diag += 4*getTwoIntIIJJ(mo1, mo2) - 2*getTwoIntIJJI(mo1, mo2);
		}
	}
	return diag;
}

EnergyType	DiagIntContainer::calcDiagClosed(const MOType *mos, INT nMOs, MOType mo1) const
{
EnergyType	e = getOneInt(mo1);
	for ( INT j=0 ; j<nMOs ; j++ )
	{
	MOType	mo2 = mos[j];
		e += 2*getTwoIntIIJJ(mo1, mo2) - getTwoIntIJJI(mo1, mo2);
	}
	return e;
}


EnergyType	DiagIntContainer::calcDiag(const Configuration<MOType> &conf) const
{
EnergyType	diag = 0;

	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
	{
	MOType	moi = conf.getOpenShell(i);
		diag += getOneInt(moi);
	}
//	cout << "diag= " << diag << endl;
	
	for ( INT k=0 ; k<conf.getNumberOfClosedShells() ; k++ )
	{
	MOType	mok = conf.getClosedShell(k);
		diag += 2*getOneInt(mok);
		diag += getTwoIntIIII(mok);
		
		for ( INT j=0 ; j<conf.getNumberOfOpenShells() ; j++ )
		{
		MOType	moj = conf.getOpenShell(j);
			diag += 2*getTwoIntIIJJ(moj, mok);
			diag -= getTwoIntIJJI(moj, mok);
		}
	}
//	cout << "diag= " << diag << endl;
	
	for ( INT j=0 ; j<conf.getNumberOfOpenShells()-1 ; j++ )
	{
	MOType	moj = conf.getOpenShell(j);
		for ( INT i=j+1 ; i<conf.getNumberOfOpenShells() ; i++ )
		{
		MOType	moi = conf.getOpenShell(i);
			diag += getTwoIntIIJJ(moi, moj);
		}
	}
//	cout << "diag= " << diag << endl;
	
	for ( INT k=0 ; k<conf.getNumberOfClosedShells()-1 ; k++ )
	{
	MOType	mok = conf.getClosedShell(k);
		for ( INT l=k+1 ; l<conf.getNumberOfClosedShells() ; l++ )
		{
		MOType	mol = conf.getClosedShell(l);
			diag += 4*getTwoIntIIJJ(mok, mol);
			diag -= 2*getTwoIntIJJI(mok, mol);
		}	
	}
//	cout << "diag= " << diag << endl;


///////////////////////////////////////////////////////////////////////
//
//	handle last part concerning representation matrices
//



	if ( conf.getNumberOfOpenShells() )
	{
	IntegralType	intBuf[MAXOPENSHELLS*(MAXOPENSHELLS+1)/2];

	INT	l = 0;

	//	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//
	//	attention: ordering of permutations from "calcdarhp5.f"!
	//
	//	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		for ( INT i=1 ; i<conf.getNumberOfOpenShells() ; i++ )
			for ( INT j=0 ; j<i ; j++ )
				intBuf[l++] = getTwoIntIJJI(conf.getOpenShell(i), conf.getOpenShell(j));



	HMatElements<MatrixType, VectorType>	*repMatP5 = repMatsP5[conf.getNumberOfOpenShells()];

	const MatrixType	*pMatP5 = repMatP5->getCbMat()->getP();
	INT	dim = repMatP5->getNumberOfRows();
	EnergyType	minE = 1e+20;
	EnergyType	sum = 0;

		for ( INT i=0 ; i<dim ; i++ )
		{
		double	h = 0;
			pMatP5 += i*l;

			for ( INT ll=0 ; ll<l ; ll++ )
				h += intBuf[ll] * *pMatP5++;

			if ( h<minE )
				minE = h;
			sum += h;
		}
//		diag += minE;
		diag += sum/dim;
	}
	

	return diag;
}



Configuration<MOType>	DiagIntContainer::calcGroundStateClosed(INT nElectrons)
{

	n = 0;
	nMOs = nElectrons/2;
	minDiag = 1e+20;

	occInd = new INT[moirreps->getNumberOfIrReps()];
	minInd = new INT[moirreps->getNumberOfIrReps()];
TimeTicks	ticks;
	cout << "searching closed shell ground state configuration..." 
		<< flush;
	ticks.start();
	dist(nMOs, moirreps->getNumberOfIrReps()-1);
	ticks.stop();
	cout << ticks << endl;
	delete occInd;
	
	minDiag += core;
	
	cout << "minDiag= " << minDiag << endl;;
	

Configuration<MOType>	conf;
	for ( IrRep i=0 ; i<moirreps->getNumberOfIrReps() ; i++ )
		for ( INT j=0 ; j<minInd[i] ; j++ )
			conf.appendClosedShell(moirreps->getStartMO(i) + j);


	delete minInd;

	return conf;
}






void	DiagIntContainer::performExcitation(
	INT toAnnihilate, INT toCreate, 
	MOType mo1, MOType mo2, 
	Configuration<MOType> conf)
{
	if ( toAnnihilate )
	{
		for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		{
		MOType	mo = conf.getOpenShell(i);

			if ( !active[mo] || mo<mo1 )
				continue;

			conf.annihilate(mo);
			performExcitation(toAnnihilate-1, toCreate, mo, 1, conf);
			conf.create(mo);
		}
		for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		{
		MOType	mo = conf.getClosedShell(i);

			if ( !active[mo] || mo<mo1 )
				continue;

			conf.annihilate(mo);
			performExcitation(toAnnihilate-1, toCreate, mo, 1, conf);
			conf.create(mo);
		}
		return;
	}
	if ( toCreate )
	{
	MOType	mos = mo2;
	MOType	moe = maxMO;
	
		for ( MOType mo=mos ; mo<=moe ; mo++ )
		{
			if ( active[mo] )
			{
				if ( !conf.create(mo) )
				{
					performExcitation(0, toCreate-1, 0, mo, conf);
					conf.annihilate(mo);
				}
			}
		}
		return;
	}

	if ( conf.getNumberOfOpenShells()<minOpenShells 
		|| conf.getNumberOfOpenShells()>maxOpenShells )
		return;

	n++;
//	cout << conf << endl;		


	TConf tc(calcDiag(conf), conf);
	lowest[conf.calcIrRep(*moirreps)].insert(tc);

}


ConfigurationSet	DiagIntContainer::calcRefConfs(
	Configuration<MOType> groundState, INT nConfs)
{
	

IrRep	irrep = groundState.calcIrRep(*moirreps);

	n = 0;
	

TimeTicks	ticks;
	cout << "searching for " << nConfs << " configurations having lowest energy" << endl;
	cout << "within single/double excitation space..." << flush;
	ticks.start();

//FD	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()](nConfs);
	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()];
	for (INT i=0; i<moirreps->getNumberOfIrReps(); i++)
	{
		lowest[i].setMaxEntries(nConfs);
	}

	for ( INT i=1 ; i<=2 ; i++ )
		performExcitation(i, i, 1, 1, groundState);
		
	ticks.stop();
	cout << ticks << endl;
	
	
	cout << n << " configurations scanned" << endl;

ConfigurationSet	set;

	for ( INT i=0 ; i<nConfs ; i++ )
	{
	TConf	tc(lowest[irrep][i]);
	Configuration<MOType>	c(tc.conf);
		set.add(c);	
		cout << tc.conf << ": " << setprecision(12) << tc.e+core << endl;
	}
	
	delete[] lowest;
	return set;
}



MOEquivalence *	DiagIntContainer::calcMODegeneration(INT nElectrons)
{
EnergyType* e = new EnergyType[maxMO];
/*
Configuration<MOType>	ground = calcGroundStateClosed(nElectrons);

	for ( INT i=0 ; i<ground.getNumberOfClosedShells() ; i++ )
	{
	Configuration<MOType>	conf(ground);
	MOType	mo = conf.getClosedShell(i);
		conf.deleteDoubleMO(mo);
				
		e[mo-1] = calcDiagClosed(conf.getClosedShellP(), nMOs-1);
		cout << mo << ", " << conf << ": " << e[mo-1] << endl;
	}
	cout << "==================" << endl;
	for ( MOType mo=1 ; mo<=maxMO ; mo++ )
	{
	INT	flag = 0;
		for ( INT i=0 ; i<ground.getNumberOfClosedShells() ; i++ )
			if ( (flag=(ground.getClosedShell(i)==mo)) )
				break;
		if ( flag )
			continue;
				
	Configuration<MOType>	conf(ground);
		conf.create(mo);
		conf.create(mo);
		
		e[mo-1] = calcDiagClosed(conf.getClosedShellP(), nMOs+1);
		cout << mo << ", " << conf << ": " << e[mo-1] << endl;
	}
*/
	for ( MOType mo=0 ; mo<maxMO ; mo++ )
	{
	Configuration<MOType>	conf;
		conf.create(mo+1);
		conf.create(mo+1);
				
		e[mo] = calcDiagClosed(conf.getClosedShellP(), 1);

//	MOType	occ[maxMO] = 
//		e[mo] = calcDiagClosed(conf.getClosedShellP(), 0, mo);

		cout << mo << ", " << conf << ": " << e[mo] << endl;
	}
	return new MOEquivalence(e, maxMO);

	delete[] e;
}





SortedList<DiagIntContainer::TConf>	*	DiagIntContainer::calcRefConfs(
		INT nElectrons,
		INT nConfs)
{
Configuration<MOType>	groundClosedShell = calcGroundStateClosed(nElectrons);

/*MOType	mo = 0;

	if ( nElectrons & 1 )
	{
		for ( INT i=0 ; i<groundClosedShell.getNumberOfClosedShells() ; i++ )
			if ( moirreps->getIrRep(groundClosedShell.getClosedShell(i))==irrep )
				mo = groundClosedShell.getClosedShell(i);

		groundClosedShell.create(mo+1);
	}
*/
//FD	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()](nConfs);
	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()];
        for (INT i=0; i<moirreps->getNumberOfIrReps(); i++)
        {
                lowest[i].setMaxEntries(nConfs);
        }

	n = 0;
	
	cout << groundClosedShell << endl;
	
	cout << endl;
	cout << "active space: " << endl;
	createActiveSpace(groundClosedShell);
	for ( INT i=1 ; i<=maxMO ; i++ )
		if ( active[i] )
			cout << i << " ";
	cout << endl;
	cout << endl;
	
TimeTicks	ticks;
	cout << "searching for " << nConfs << " configurations having lowest energy" << endl;
	cout << "within single/double/triple/quadruple excitation space..." << flush;
	ticks.start();

//	performExcitation(1, 1, 1, 1, groundClosedShell);
//	for ( INT i=1 ; i<=4 ; i++ )
	for ( INT i=1 ; i<=2 ; i++ )
		performExcitation(i - (nElectrons & 1), i, 1, 1, groundClosedShell);
	
	ticks.stop();
	cout << ticks << endl;
	
	
	cout << n << " configurations scanned" << endl;
	
/*ConfigurationSet	set;

	for ( INT i=0 ; i<nConfs ; i++ )
		cout << (*lowest)[i].conf << ": " << setprecision(12) << (*lowest)[i].e+core << endl;

	if ( killDegen )
	{
		cout << endl << "using: " << endl;
		for ( INT i=0 ; i<nConfs ; i++ )
			if ( fabs((*lowest)[i].e-(*lowest)[0].e)<1e-6 )
			{
			TConf	tc((*lowest)[i]);
			Configuration<MOType>	c(tc.conf);
				set.add(c);	
				cout << tc.conf << ": " << setprecision(12) << tc.e+core << endl;
			}
		cout << endl << endl;
	}
	else
	{
		for ( INT i=0 ; i<nConfs ; i++ )
		{
		Configuration<MOType>	c((*lowest)[i].conf);
			set.add(c);	
		}
	}
*/	
	return lowest;
}


ConfigurationSet	DiagIntContainer::calcRefConfs(
		INT nElectrons,
		IrRep _irrep,
		Configuration<MOType> &groundState,
		INT nConfs,
		INT killDegen)
{
Configuration<MOType>	groundClosedShell = calcGroundStateClosed(nElectrons);


	irrep = _irrep;
/*MOType	mo = moirreps->getStartMO(irrep) - 1;


	if ( nElectrons & 1 )
	{
		for ( INT i=0 ; i<groundClosedShell.getNumberOfClosedShells() ; i++ )
			if ( moirreps->getIrRep(groundClosedShell.getClosedShell(i))==irrep )
				mo = groundClosedShell.getClosedShell(i);

		groundClosedShell.create(mo+1);
	}
*/

	
	n = 0;
	
//FD	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()](nConfs);
	lowest = new SortedList<TConf>[moirreps->getNumberOfIrReps()];
        for (INT i=0; i<moirreps->getNumberOfIrReps(); i++)
        {
                lowest[i].setMaxEntries(nConfs);
        }
	
	cout << groundClosedShell << endl;
	
	cout << endl;
	cout << "active space: " << endl;
	createActiveSpace(groundClosedShell);
	for ( INT i=1 ; i<=maxMO ; i++ )
		if ( active[i] )
			cout << i << " ";
	cout << endl;
	cout << endl;
	
TimeTicks	ticks;
	cout << "searching for " << nConfs << " configurations having lowest energy" << endl;
	cout << "within single/double excitation space..." << flush;
	ticks.start();

//	performExcitation(1, 1, 1, 1, groundClosedShell);
	for ( INT i=1 ; i<=2 ; i++ )
//	for ( INT i=1 ; i<=4 ; i++ )
		performExcitation(i - (nElectrons & 1), i, 1, 1, groundClosedShell);
	
	ticks.stop();
	cout << ticks << endl;
	
	
	cout << n << " configurations scanned" << endl;
	
ConfigurationSet	set;


TConf* confs = new TConf[nConfs];

	{
	INT	j = 0;
		for ( SortedList<DiagIntContainer::TConf>::const_iterator 
			i=lowest[irrep].begin() ; i!=lowest[irrep].end() ; ++i )
		{
			confs[j++] = *i;
		}
	}

	groundState = confs[0].conf;

	for ( INT i=0 ; i<nConfs ; i++ )
		cout << confs[i].conf << ": " << setprecision(12) << confs[i].e+core << endl;

	if ( killDegen )
	{
		cout << endl << "using: " << endl;
		for ( INT i=0 ; i<nConfs ; i++ )
			if ( fabs(confs[i].e-confs[0].e)<1e-6 )
			{
			Configuration<MOType>	c(confs[i].conf);
				set.add(c);	
				cout << confs[i].conf << ": " << setprecision(12) << confs[i].e+core << endl;
			}
		cout << endl << endl;
	}
	else
	{
		for ( INT i=0 ; i<nConfs ; i++ )
		{
		Configuration<MOType>	c(confs[i].conf);
			set.add(c);	
		}
	}

	delete[] confs;
	delete[] lowest;	
	return set;

}

void	DiagIntContainer::createActiveSpace(Configuration<MOType> conf)
{
INT* highest = new INT[nIrreps];
bool* occ = new bool[nIrreps];

	for ( IrRep i=0 ; i<nIrreps ; i++ )
	{
		highest[i] = moirreps->getStartMO(i) - 1;
		occ[i] = false;
	}
	
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
	{
		highest[moirreps->getIrRep(conf.getOpenShell(i))] =
			conf.getOpenShell(i);
		occ[moirreps->getIrRep(conf.getOpenShell(i))] = true;
	}
	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
	{
		if ( highest[moirreps->getIrRep(conf.getClosedShell(i))] < 
					conf.getClosedShell(i) )
			highest[moirreps->getIrRep(conf.getClosedShell(i))] =
				conf.getClosedShell(i);
		occ[moirreps->getIrRep(conf.getClosedShell(i))] = true;
	}
				
	for ( IrRep i=0 ; i<nIrreps ; i++ )
	{
//		if ( occ[i] )
		{
			for ( MOType j=highest[i] ; 
				j>highest[i]-activeOccupied && j>0 && moirreps->getIrRep(j)==i ; j-- )
				active[j] = 1;

			for ( MOType j = highest[i]+1 ;
				j<=highest[i]+activeVirtual && j<=maxMO && moirreps->getIrRep(j)==i ; j++ )
				active[j] = 1;
		}
	}

delete highest;
delete occ;
				
}




void	DiagIntContainer::dist(MOType nMOsInIrRep, IrRep nSym)
{
	if ( nMOsInIrRep<0 )
	{
		return;
	}


//	cout << " " << setw(nSym*4) << nMOsInIrRep << setw(4) << endl;
	if ( !nSym ) 
	{
		n++;
		occInd[0] = nMOsInIrRep;

	EnergyType	diag = calcDiagClosed(occInd);

		if ( diag<minDiag )
		{
			minDiag = diag;
			memcpy(minInd, occInd, nIrreps*sizeof(INT));
		}

		return;
	}

INT	h = nMOsInIrRep-moirreps->getInIrRep(nSym);
	if ( h<0 )
		h = 0;
	for ( INT i=nMOsInIrRep-h ; i>=0 ; i-- )
	{
		occInd[nSym] = i;
		dist(nMOsInIrRep-i, nSym-1);
	}
}

