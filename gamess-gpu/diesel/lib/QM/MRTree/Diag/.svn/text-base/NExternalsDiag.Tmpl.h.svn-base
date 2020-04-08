//***********************************************************************
//
//	Name:			NExternalsDiag.Tmpl.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			07.10.1998
//
//
//
//
//
//***********************************************************************

#ifndef __NExternalsDiag_Tmpl_H
#define __NExternalsDiag_Tmpl_H

#include "../../../../config.h"

#include "NExternalsDiag.h"
#include "InternalConfsDiag.h"
#include "TupelStructureDiag.h"
#include "extMOsDiag.h"

#include <errno.h>

#include "../../IO/TimeTicks.h"
#include "../../IO/Verbosity.h"
#include "../../MRCIMatrix/InteractionClasses.h"
#include "MatchingExternalMOIterator.h"
#include "../../IntegralIndex/TwoElectronIntegralIndexUpdateMask.h"
#include "../../RepresentationMatrices/RepresentationMatrixFortranInterface.h"
#include "../../Cache/VarSizeReadOnlyCache.h"
#include "../../Cache/TableKey.h"
#include "../MatrixAction.h"
#include "../../../Math/etc/BinomialCoefficient.h"


#include "../../MRCIMatrix/MatrixAction/ActionMultiplication.h"
#include "../../MRCIMatrix/MatrixAction/ActionDensity.h"
#include "../../MRCIMatrix/MatrixAction/ActionMatrixStorage.h"
#include "../../MRCIMatrix/MatrixAction/ActionPointerStorage.h"

#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>

#ifdef _AIX
#include <strings.h>
#include <sys/select.h>
#endif

using std::flush;

template <class MatrixType, class VectorType, class MatrixAction, class MatrixElementCalculator>
class TwoDimIterator {
public:
	TwoDimIterator(
		INT numberOfSlaves,
		const	MRMOs	*mrmos,
		SharedMemory	*Int4sharedMem,			// shared memory containing 4 index integrals
		const NExternalsDiag *mrccA,
		const NExternalsDiag *mrccB,
		MatrixAction &action,
		const MatrixElementCalculator &dummy,
		INT highestExcitation,
		INT	upperTriangular);

	~TwoDimIterator();

	
	CIVectors<VectorType> *	getXVector()
	{	return xShared;	}

	CIVectors<VectorType> *	getYVector()
	{	return yShared[0];	}

//	VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	getCache() const
//	{	return cache;	}


struct TJob {
	TJob() {}
	
	TJob(INT _areaA, INT _areaB, INT _indStart, INT _indEnd)
	{
		areaA = _areaA;
		areaB = _areaB;
		indStart = _indStart;
		indEnd = _indEnd;
	}
	
	INT	areaA;
	INT	areaB;
	INT	indStart;
	INT	indEnd;
	};


private:
	void	innerIteration(
		INT	ia,
		INT	ib,
		INT	ja,
		INT	jb,
		const TupelStructureDiag	*mainA,
		const TupelStructureDiag	*mainB);
	void	receiveOrders();


DiffConf<MOType>	diffConfExtAInternalB, diffConfInternal;
BinomialCoefficient	*binom;
InteractionClasses<MatrixElementCalculator>	*classes;
VarSizeReadOnlyCache<TableKey, MatrixElementCalculator>	*cache;

INT	numberOfSlaves;
const	MRMOs	*mrmos;
SharedMemory	*Int4sharedMem;			// shared memory containing 4 index integrals
const NExternalsDiag *mrccA;
const NExternalsDiag *mrccB;
MatrixAction &action;
INT highestExcitation;
INT	upperTriangular;

//------------------------------------------------------------------------
//
//		inter process communication
//
//------------------------------------------------------------------------


INT	isSlave;							// flag if instance belongs to
										// slave process

CIVectors<VectorType>	*xShared;		// CI-vectors in shared memory
SharedMemory	*ySharedMem;			// shared memory containing CI-y vectors
CIVectors<VectorType>	**yShared;		// CI-vectors in shared memory


struct TPipe {
	int out;
	int in;
	};
struct TComm {
	TPipe	toSlave;
	TPipe	toMaster;
	};
	
TPipe	*comm;							// communication pipes to slaves

};



//********************************************************************
// bug in GNU C/C++ 2.8.1:
// the following macors generate an internal compiler error
// when called within MRCIMatrix<MatrixType, VectorType>

static void	Local_fd_zero(fd_set *set)
{
	FD_ZERO(set);
}

static void	Local_fd_set(INT i, fd_set *set)
{
	FD_SET(i, set);
}

static void	Local_fd_clr(INT i, fd_set *set)
{
	FD_CLR(i, set);
}

static INT	Local_fd_isset(INT i, fd_set *set)
{
	return FD_ISSET(i, set);
}

//
//
//********************************************************************

template <class MatrixType, class VectorType, class MatrixAction, class MatrixElementCalculator>
TwoDimIterator<MatrixType, VectorType, MatrixAction, MatrixElementCalculator>::TwoDimIterator(
	INT numberOfSlaves,
	const	MRMOs	*mrmos,
	SharedMemory	*Int4sharedMem,
	const NExternalsDiag *mrccA,
	const NExternalsDiag *mrccB,
	MatrixAction &action,
	const MatrixElementCalculator &dummy,
	INT highestExcitation,
	INT	upperTriangular) :
		numberOfSlaves(numberOfSlaves),
		mrmos(mrmos),
		Int4sharedMem(Int4sharedMem),
		mrccA(mrccA),
		mrccB(mrccB),
		action(action),
		highestExcitation(highestExcitation),
		upperTriangular(upperTriangular)
{
TimeTicks	ticks, ticksTotal;
TimeTicks	ticksSum;

	
	binom = new BinomialCoefficient(MAXOPENSHELLS);

	classes = new InteractionClasses<MatrixElementCalculator>(binom);

	cache = new VarSizeReadOnlyCache<TableKey, MatrixElementCalculator >
			(1024, (1<<20));


//**************************************************************************
//**************************************************************************
//**************************************************************************



	xShared = NULL;
	yShared = NULL;
	ySharedMem = NULL;



INT	roots = 0;
		roots = mrccA->getNumberOfRoots();
INT	totalDim = mrccA->getNumberOfTotalSpinAdaptedFunctions();
	if ( !roots || !totalDim )
		numberOfSlaves = 0;
	if ( numberOfSlaves )
	{
//		roots = ((ActionMultiplication<MatrixType, VectorType> &) action).CIy.getN();
//		cout << "roots=" << roots << endl;
//		cout << "totalDim=" << totalDim << endl;
		xShared = new CIVectors<VectorType>(totalDim, roots, SharedMemory::Master);
		
		memcpy(xShared->getP(), ((ActionMultiplication<MatrixType, VectorType> &) action).CIx.getP(),
			totalDim*roots*sizeof(VectorType));
		
		ySharedMem = new SharedMemory(
			numberOfSlaves*totalDim*roots*sizeof(VectorType), SharedMemory::Master);
		yShared = new CIVectors<VectorType> * [numberOfSlaves];
		for ( INT i=0 ; i<numberOfSlaves ; i++ )
			yShared[i] = new CIVectors<VectorType>(
				totalDim, roots, (VectorType *) ySharedMem->allocate(totalDim*sizeof(VectorType), roots));
	}


struct TComm Comm;
	


	isSlave = 0;
	
	comm = new TPipe[numberOfSlaves];
	
	for ( INT i=0 ; i<numberOfSlaves ; i++ )
	{
		pipe((int *) &Comm.toSlave);
		pipe((int *) &Comm.toMaster);
	
//		cout << "forking child process #" << i+1 << "..." << flush;
		
	pid_t	pid = fork();
		if ( pid==-1 )
		{
			cout << "fork failed" << endl;
			exit(1);
		}
//		cout << endl;
		if ( pid )
		{	// parent process
			comm[i].in = Comm.toMaster.out;
			close(Comm.toMaster.in);
			comm[i].out = Comm.toSlave.in;
			close(Comm.toSlave.out);
		}
		else
		{	// child process
			isSlave = i+1;
//			cout << "I am child #" << isSlave << endl;
			comm->in = Comm.toSlave.out;
			close(Comm.toSlave.in);
			comm->out = Comm.toMaster.in;
			close(Comm.toMaster.out);
			
			
			if ( Int4sharedMem )
				Int4sharedMem->setSlave();
			if ( xShared )
				xShared->setSlave();

//			xv = new CIVectors<VectorType>(totalDim, 1, xShared->getP());
//			yv = new CIVectors<VectorType>(totalDim, roots, yShared[isSlave]->getP());

			delete ((ActionMultiplication<MatrixType, VectorType> &) action).CIx.getP();
			delete ((ActionMultiplication<MatrixType, VectorType> &) action).CIy.getP();
//			((ActionMultiplication<MatrixType, VectorType> &) action).CIx.setP(xShared->getP());
//			((ActionMultiplication<MatrixType, VectorType> &) action).CIy.setP(yShared[isSlave-1]->getP());
			((ActionMultiplication<MatrixType, VectorType> &) action).CIx = *xShared;
			((ActionMultiplication<MatrixType, VectorType> &) action).CIy = *yShared[isSlave-1];


//			action.CIy = *yv;

//				ActionMultiplication<MatrixType, VectorType> a(*xv, *yv);
//			((ActionMultiplication<MatrixType, VectorType> &) action) = a;
			break;
		}
	}
	
	
	if ( isSlave )
	{
		receiveOrders();
//		delete this;
//		cout << "child #" << isSlave << " dies." << endl;
		_exit(0);
	}


//**************************************************************************
//**************************************************************************
//**************************************************************************




fd_set	test_set, ready_set;
INT	max_fd = numberOfSlaves ? comm[numberOfSlaves-1].in+1 : 0;

INT* out_fd = new INT[max_fd];
INT* ready = new INT[max_fd];
char	buf[1000];

	for ( INT i=0 ; i<max_fd ; i++ )
		ready[i] = 0;

	Local_fd_zero(&test_set);
	for ( INT i=0 ; i<numberOfSlaves ; i++ )
	{
		Local_fd_set(comm[i].in, &test_set);
		out_fd[comm[i].in] = comm[i].out;
		ready[comm[i].in] = 0;
//		cout << "slave: " << i << " <--> pipe: " << comm[i].in << endl;
	}


	for ( INT i=0 ; i<numberOfSlaves ; i++ )
		yShared[i]->clear();






	ticksTotal.start();


	mrmos = (MRMOs *) mrccA->getMRMOs();


//	main loops over intern-ia / intern-ib
	for ( INT ia=0 ; ia<mrccA->getNumberOfElements() ; ia++ )
	{
	const InternalConfsDiag	*mainContA = (*mrccA)[ia];

//	calculate upper triangular matrix only
		for ( INT ib=(upperTriangular ? ia : 0) ; ib<mrccB->getNumberOfElements() ; ib++ )
		{
//	no interaction if difference in more than 2 MOs
			if ( abs(ia-ib)>highestExcitation )
				continue;

			if ( verbosity.isActive(Verbosity::IterationBlocks) )
			{
				cout << "<intern-" << ia << "|H|intern-" << ib << ">...";
				cout.flush();
			}
			ticksSum.clear();
		const InternalConfsDiag	*mainContB = (*mrccB)[ib];




//*******************************************************************************
//*******************************************************************************
//*******************************************************************************


			if ( !numberOfSlaves )
			{
				for ( INT ja=0 ; ja<mainContA->getNumberOfElements() ; ja++ )
				{
					ticks.start();
				const TupelStructureDiag	*mainA = (*mainContA)[ja];		

	//	if not density calculation generate upper triangular matrix only
					for ( INT jb=(upperTriangular ? (ia==ib ? ja : 0 ) : 0) ; 
						jb<mainContB->getNumberOfElements() ; jb++ )
					{

					const TupelStructureDiag	*mainB = (*mainContB)[jb];		


						innerIteration(ia, ib, ja, jb, mainA, mainB);

					} // internal-conf B
					ticks.stop();
					ticksSum += ticks;
				} // internal-conf A
			}
			else
			{
			const INT granularity = 30;
			
				ticksSum.clear();

			INT n = (ia==ib && upperTriangular) ? 
				mainContA->getNumberOfElements()*(mainContA->getNumberOfElements()+1)/2 :
				mainContA->getNumberOfElements()*mainContB->getNumberOfElements();

			INT	step = n / granularity + 1;
				for ( INT j=0 ; j<n ; j+=step )
				{
					ticks.start();

				TJob	job(ia, ib, j, (j+step-1<n) ? j+step-1 : n-1);
					while ( 1 )
					{
					INT	breaked = 0;
						for ( INT fd=0 ; fd<max_fd ; fd++ )
							if ( ready[fd] )
							{
								read(fd, buf, 1000);
								write(out_fd[fd], &job, sizeof(TJob));
								ready[fd] = 0;
								breaked = 1;
								break;
							}
						if ( breaked )
							break;

						memcpy(&ready_set, &test_set, sizeof(test_set));
						select(max_fd+1, ((fd_set *) &ready_set), NULL, NULL, NULL);
						for ( INT fd=0 ; fd<=max_fd ; fd++ )
							if ( Local_fd_isset(fd, &ready_set) )				
								ready[fd] = 1;
					}
					ticks.stop();
					ticksSum += ticks;
				}
			}


//*******************************************************************************
//*******************************************************************************
//*******************************************************************************


			if ( verbosity.isActive(Verbosity::IterationBlocks) )
				cout << ticksSum << endl;
		} // internal-n B
	} // internal-n A
	ticksTotal.stop();
	if ( verbosity.isActive(Verbosity::IterationBlocks) )
		cout << "total: " << ticksTotal << endl;


	if ( numberOfSlaves )
	{	
		if ( verbosity.isActive(Verbosity::IterationBlocks) )
			cout << "waiting for slaves to finish..." << flush;
	//	wait for all slaves to finish
		
		for ( INT i=0 ; i<numberOfSlaves ; )
		{
		INT	fd = comm[i].in;
			if ( ready[fd] )
			{
				i++;
				continue;
			}
				
//			cout << "master: waiting for slave #" << i+1 << "..." << flush; 
				

			memcpy(&ready_set, &test_set, sizeof(test_set));
			select(max_fd+1, ((fd_set *) &ready_set), NULL, NULL, NULL);
			for ( INT ifd=0 ; ifd<=max_fd ; ifd++ )
				if ( Local_fd_isset(ifd, &ready_set) )				
				{
					ready[ifd] = 1;
					Local_fd_clr(ifd, &test_set);
				}
					
//			cout << "finished." << endl;
		}
		cout << "finished." << endl;


		if ( verbosity.isActive(Verbosity::IterationBlocks) )
			cout << "collecting results..." << flush;

	INT	roots = ((ActionMultiplication<MatrixType, VectorType> &) action).CIy.getN();
	//	collect results
		for ( INT i=1 ; i<numberOfSlaves ; i++ )
		{
		VectorType	*p0 = yShared[0]->getP();
		VectorType	*pi = yShared[i]->getP();
			for ( INT j=0 ; j<totalDim*roots ; j++ )
				*p0++ += *pi++;
		}
/*		for ( INT j=0 ; j<totalDim*roots ; j++ )
			cout << yShared[0]->getP()[j] << "\t\t" <<
			((ActionMultiplication<MatrixType, VectorType> &) action).CIy.getP()[j] <<
				endl;
*/

		memcpy( ((ActionMultiplication<MatrixType, VectorType> &) action).CIy.getP(), yShared[0]->getP(),
				totalDim*roots*sizeof(VectorType));

		if ( verbosity.isActive(Verbosity::IterationBlocks) )
			cout << "finished." << endl;
	}
	delete[] ready;
	delete[] out_fd;
}


template <class MatrixType, class VectorType, class MatrixAction, class MatrixElementCalculator>
TwoDimIterator<MatrixType, VectorType, MatrixAction, MatrixElementCalculator>::~TwoDimIterator()
{
	delete cache;
	delete classes;
	delete binom;
	if ( numberOfSlaves )
	{

		if ( !isSlave )
		{
			if ( xShared )
				delete xShared;
				
			for ( INT i=0 ; i<numberOfSlaves ; i++ )
			{
				close(comm[i].in);
				close(comm[i].out);
			}

			for ( INT i=0 ; i<numberOfSlaves ; i++ )
				wait(0);

			if ( ySharedMem )
				delete ySharedMem;
			if ( yShared )
				delete yShared;
		}
		else
		{
			close(comm->in);
			close(comm->out);
			if ( yShared )
				if ( yShared[isSlave-1] )
					delete yShared[isSlave-1];
		}
		delete comm;
	}
}





template <class MatrixType, class VectorType, class MatrixAction, class MatrixElementCalculator>
void	TwoDimIterator<MatrixType, VectorType, MatrixAction, MatrixElementCalculator>::receiveOrders()
{
TJob	job;

	// send ready signal
	write(comm->out, &job, 1);	
	while ( read(comm->in, &job, sizeof(TJob))>0 )
	{
//		cout << "child #" << isSlave << " received message from parent:" << endl;
//		cout << job.areaA << " " << job.areaB << " " 
//			<< job.indStart << " " << job.indEnd  << endl;


	const InternalConfsDiag	*mainContA = (*mrccA)[job.areaA];
	const InternalConfsDiag	*mainContB = (*mrccA)[job.areaB];

		for ( INT i=job.indStart ; i<=job.indEnd ; i++ )
		{
		INT	noA, noB;

			if ( job.areaA==job.areaB && upperTriangular )
			{
				noB = (INT) floor(-0.5 + sqrt(0.25+2*i));
				noA = i - noB*(noB+1)/2;
			}
			else
			{
				noA = i / mainContB->getNumberOfElements();
				noB = i % mainContB->getNumberOfElements();
			}

//			cout << "::" << job.areaA << " " << job.areaB << " " 
//				<< noA << " " << noB << endl;


		const TupelStructureDiag	*mainA = (*mainContA)[noA];
		const TupelStructureDiag	*mainB = (*mainContB)[noB];

			innerIteration(job.areaA, job.areaB, noA, noB, mainA, mainB);

		}


		// send ready signal
		write(comm->out, &job, 1);	
	}
}






template <class MatrixType, class VectorType, class MatrixAction, class MatrixElementCalculator>
void	TwoDimIterator<MatrixType, VectorType, MatrixAction, MatrixElementCalculator>::
	innerIteration(
	INT	ia,
	INT	ib,
	INT	ja,
	INT	jb,
	const TupelStructureDiag	*mainA,
	const TupelStructureDiag	*mainB)
{

										

//		no interaction if more than highestExcitation excitation between the two main-n
INT	orderInt = 0;

	if ( ib>=ia )
	{
		if ( (orderInt=Configuration<MOType>::calcExcitationOrderFast(
			*mainA, *mainB))>highestExcitation )
			return;
	}
	else
	{
		if ( (orderInt=Configuration<MOType>::calcExcitationOrderFast(
			*mainB, *mainA))>highestExcitation )
			return;
	}

INT	maxExcitation = 0;
	if ( ib>ia )
		maxExcitation = orderInt - (ib-ia) + ib;
	else
		maxExcitation = orderInt - (ia-ib) + ia;


	diffConfInternal.calcDiffConf(*mainA, *mainB);


	for ( INT ka=0 ; ka<mainA->getNumberOfElements() ; ka++ )
	{
	MOType	extA[2] = { 0, 0 };
	const extMOsDiag	*MOA = (*mainA)[ka];
		if ( !MOA )
			continue;

	INT	safsA = MOA->getSAFInc();





		for ( INT la=0 ; la<MOA->getNumberOfElements() ; la++ )
		{
		INT safStartA = MOA->getSAFStart() + la*safsA;

			for ( INT ii=0 ; ii<MOA->getNumberOfTotalMOs() ; ii++ )
				extA[ii] = (*MOA)[la][ii];

//---------------------------------------------------
			diffConfExtAInternalB = diffConfInternal;

			if ( MOA->getNumberOfTotalMOs() )
				diffConfExtAInternalB.addExternal(
					Configuration<MOType>(
						MOA->getNumberOfOpenMOs(),
						MOA->getNumberOfClosedMOs(),
						(*MOA)[la]),
						Configuration<MOType>());

		MOType	extB[2] = { 0, 0 };




	INT	jDiag = ib==ia && jb==ja && upperTriangular;
			for ( INT kb=(jDiag) ? ka : 0 ; 
				kb<mainB->getNumberOfElements() ; kb++ )
			{
			const extMOsDiag	*MOB = (*mainB)[kb];
				if ( !MOB )
					continue;
			INT	safsB = MOB->getSAFInc();
			INT	openMOs = MOB->getNumberOfOpenMOs();
			INT	OldCmpStatus = -1;
			INT	QOldCmpStatus = -1;
			INT	kDiag = jDiag && kb==ka;
			IrRep	oldIrrep = 1000;
			typename InteractionClasses<MatrixElementCalculator>::TCase	*actualClass = NULL;


			if ( !MOB->getNumberOfElements() )
				continue;

			classes->clear();

			MatchingExternalMOIterator iter(MOA, MOB, 
					maxExcitation, highestExcitation, kDiag, la);
				while ( iter.next() ) 
				{
				INT	lb = iter.getIndex();


					for ( INT ii=0 ; ii<MOB->getNumberOfTotalMOs() ; ii++ )
						extB[ii] = (*MOB)[lb][ii];

					if ( MOB->getNumberOfTotalMOs() && mrmos->getIrRep(extB[0])!=oldIrrep )
					{
						oldIrrep = mrmos->getIrRep(extB[0]);
						classes->clear();
					}


				INT	CmpStatus = 
						((extA[0]==extB[0]) << 0) |
						((extA[1]==extB[1]) << 1) |
						((extA[0]==extB[1]) << 2) |
						((extA[1]==extB[0]) << 3);

				INT	QCmpStatus = 
						((extA[0] <extB[0]) << 0) |
						((extA[1] <extB[1]) << 1) |
						((extA[0] <extB[1]) << 2) |
						((extA[1] <extB[0]) << 3);

					if ( CmpStatus!=OldCmpStatus )
						actualClass = &(*classes)[CmpStatus];

					// ===================
					// new total case
					// ===================
					if ( actualClass->repMats == (MatrixElementCalculator *) 1 )
						continue;


					if ( !actualClass->repMats )
					{


recalc:									actualClass->diffConf = diffConfExtAInternalB;

						actualClass->clearExtPosB();

						if ( MOB->getNumberOfTotalMOs() )
						{
						INT	orderExt = actualClass->diffConf.addExternalTo(
								openMOs, MOB->getNumberOfClosedMOs(), (*MOB)[lb]);




							for ( INT ii=0 ; ii<actualClass->to->getNumberOfOpenShells() ; ii++ )
							{
								if ( extB[0] == actualClass->to->getOpenShell(ii) )
									actualClass->extPosBOpen[0] = ii;
								if ( extB[1] == actualClass->to->getOpenShell(ii) )
									actualClass->extPosBOpen[1] = ii;
							}

							for ( INT ii=0 ; ii<actualClass->to->getNumberOfClosedShells() ; ii++ )
							{
								if ( extB[0] == actualClass->to->getClosedShell(ii) )
									actualClass->extPosBClosed[0] = ii;
								if ( extB[1] == actualClass->to->getClosedShell(ii) )
									actualClass->extPosBClosed[1] = ii;
							}


							orderExt -= ib-ia;




							if ( orderInt+orderExt > 2 )
							{
								actualClass->repMats = (MatrixElementCalculator *) 1;
								continue;
							}

							actualClass->tablecase->calcLess3(actualClass->diffConf);
							actualClass->updateMask = TwoElectronIntegralIndexUpdateMask(
								actualClass->tablecase->getCbExIndex(),
								(*MOB)[lb], MOB->getNumberOfTotalMOs());
						}
						else
						{
							actualClass->tablecase->calcLess3(actualClass->diffConf);
							actualClass->updateMask = TwoElectronIntegralIndexUpdateMask(
								actualClass->tablecase->getCbExIndex(),
								NULL, 0);
						}

						if ( actualClass->tablecase->getP()==0 )
						{
							actualClass->repMats = (MatrixElementCalculator *) 1;
							continue;
						}

						actualClass->repMats = (MatrixElementCalculator *) (*cache)[TableKey(*actualClass->tablecase)];
					}
					// ============================================================================


					else
					{
						// ====================
						// q-case changed only
						// ====================
						if ( QCmpStatus!=QOldCmpStatus || CmpStatus!=OldCmpStatus )
						{
							actualClass->diffConf.updatePos(*(actualClass->to), (INT *) actualClass->diffConf.getPosToP());
							actualClass->tablecase->updateqR(actualClass->diffConf);
							actualClass->repMats = (MatrixElementCalculator *) (*cache)[TableKey(*actualClass->tablecase)];

						}
						// ======================================================================



						//	update MOs in tablecase
						actualClass->updateMask.useOn(actualClass->tablecase->getCbExIndex(), (*MOB)[lb]);
					TwoElectronIntegralCbExIndex CbExIndex(actualClass->tablecase->getCbExIndex());

						//	attention: ordering may be wrong ==> recalculate
						if ( CbExIndex[0]<CbExIndex[1] ||
							CbExIndex[1]<CbExIndex[2] ||
							CbExIndex[2]<CbExIndex[3] )
							goto recalc;


						//	update MOs in DiffConf
						if ( actualClass->extPosBOpen[0]>=0 )
							actualClass->to->setOpenShell(actualClass->extPosBOpen[0], extB[0]);
						if ( actualClass->extPosBOpen[1]>=0 )
							actualClass->to->setOpenShell(actualClass->extPosBOpen[1], extB[1]);

						if ( actualClass->extPosBClosed[0]>=0 )
							actualClass->to->setClosedShell(actualClass->extPosBClosed[0], extB[0]);
						if ( actualClass->extPosBClosed[1]>=0 )
							actualClass->to->setClosedShell(actualClass->extPosBClosed[1], extB[1]);
					}

					OldCmpStatus = CmpStatus;
					QOldCmpStatus = QCmpStatus;

	INT safStartB = MOB->getSAFStart() + lb*safsB;


//==============================================================================================

					action.doit(
						*actualClass->repMats,
						safStartA, safStartB, 
						actualClass->diffConf, *actualClass->tablecase);


//==============================================================================================


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
				} // external MOs B
//								if ( countOK!=countTest )
//									cout << "---------- " << countOK << " " << countTest << endl;
			} // tupel structure B
		} // external MOs A
	} // tupel structure A
}





#endif
