//***********************************************************************
//
//	Name:			MRCIMatrix.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.06.1996
//					21.02.1998
//
//
//
//
//
//***********************************************************************




#include "MRCIMatrix.h"


#include "InteractionClasses.h"


#include "../MRTree/Diag/NExternalsDiag.h"
#include "../MRTree/Diag/InternalConfsDiag.h"
#include "../MRTree/Diag/TupelStructureDiag.h"
#include "../MRTree/Diag/extMOsDiag.h"
#include "../MRTree/Diag/MatchingExternalMOIterator.h"

#include "../IntegralIndex/TwoElectronIntegralIndexUpdateMask.h"

#include "../RepresentationMatrices/RepresentationMatrixFortranInterface.h"
#include "../Cache/VarSizeReadOnlyCache.h"

#include "../MRTree/Diag/NExternalsDiag.Tmpl.h"

#include "../IO/TimeTicks.h"


#include "../MRTree/MRTreeIterator.h"
#include "../Parallel/SharedMemory.h"
#include "../../Container/DiskBuffer.h"
#include "../IO/Verbosity.h"






template <class MatrixType, class VectorType>
MRCIMatrix<MatrixType, VectorType>::MRCIMatrix(NExternalsDiag	* _mrcc,
		Fort31File _Fort31,
		INT cacheEntries, INT cacheMemory, INT _numberOfSlaves) :
		CICalculation<MatrixType, VectorType>(_mrcc->getMultiplicity(),
		_numberOfSlaves, *_mrcc->getMRMOs(), _Fort31, this->storeAll),
		DavidsonMatrix<MatrixType, VectorType>(_mrcc->getNumberOfRefConfSpinAdaptedFunctions(),
		_mrcc->getNumberOfTotalSpinAdaptedFunctions())
		
{
	mrccA = mrccB = _mrcc;
	numberOfSlaves = _numberOfSlaves;
			
	GeneralizedMO::mrmos = (MRMOs *) mrccA->getMRMOs();

	cache = new VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >
			(cacheEntries, cacheMemory);

//--------------------------------------------------------------------------

	// initialize RefSAFIndex
	{
		RefSAFIndex = new INT[this->refDim];

	InternalConfsDiag *internal = mrccA->operator [] (0);
	ContainerIterator iteri = internal->first();
	INT	SAFNrA = 0;
	INT	ind = 0;
		for ( INT i=0 ; i<mrccA->getNumberOfReferences() ; i++ )
		{
			while ( !(*internal)[iteri]->isReference() )
				internal->next(iteri);

			ConfigurationSAFNr<MOType> mainA =
				(*(*internal)[iteri])[0]->getConfigurationSAFNr(0);
		INT safsA = mainA.getSAFInc();
			for ( INT k=0 ; k<safsA ; k++ )
				RefSAFIndex[ind++] = mainA.getSAFNr() + k;

			internal->next(iteri);
			SAFNrA += safsA;
		}
	}
}



template <class MatrixType, class VectorType>
MRCIMatrix<MatrixType, VectorType>::MRCIMatrix(
		NExternalsDiag *_mrccA, 
		NExternalsDiag *_mrccB,
		Fort31File _f31,
		INT cacheEntries, INT cacheMemory, INT _numberOfSlaves):
		CICalculation<MatrixType, VectorType>(_mrccA->getMultiplicity(),
		_numberOfSlaves, *_mrccA->getMRMOs(), _f31, this->storeNone),
		DavidsonMatrix<MatrixType, VectorType>(0, 0)
{
	mrccA = _mrccA;
	mrccB = _mrccB;
	numberOfSlaves = _numberOfSlaves;
			
	GeneralizedMO::mrmos = (MRMOs *) mrccA->getMRMOs();

	cache = new VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >
			(cacheEntries, cacheMemory);

	this->refDim = this->totalDim = 0;

	RefSAFIndex = NULL;
}



template <class MatrixType, class VectorType>
MRCIMatrix<MatrixType, VectorType>::MRCIMatrix(
	const CICalculation<MatrixType, VectorType> &ciCalc, NExternalsDiag *_mrcc) :
	CICalculation<MatrixType, VectorType>(ciCalc),
		DavidsonMatrix<MatrixType, VectorType>(0, 0)
{
	mrccA = mrccB = _mrcc;

	numberOfSlaves = 0;

	cache = new VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >
			(2048, (1<<20));
			
	this->refDim = this->totalDim = 0;
//	GeneralizedMO::mrmos = (MRMOs *) mrccA->getMRMOs();
	RefSAFIndex = NULL;

	this->freeOnDestruct = 0;
}



/*
MRCIMatrix<double, double>::MRCIMatrix(
		NExternalsDiag *_mrcc, 
		NExternalsDiag *_mrcc2,
		Fort31File _f31,
		INT cacheEntries, INT cacheMemory, INT _numberOfSlaves) :
		CICalculation(mrcc->getMultiplicity(), cacheEntries, cacheMemory,
		_numberOfSlaves, _mrcc->getMRMOs(), _f31)
{
	mrcc = _mrcc;
	mrcc2 = _mrcc2;
	Fort31 = _f31;
	numberOfSlaves = _numberOfSlaves;
			
	intContainer = NULL;
	this->Int4sharedMem = NULL;
	
	this->refDim = this->totalDim = 0;

	GeneralizedMO::mrmos = (MRMOs *) mrcc->getMRMOs();

	repMatFInt = 
		new RepresentationMatrixFortranInterface<MatrixType>(mrcc->getMultiplicity());

RepresentationMatrices<MatrixType, VectorType>  init1(mrcc->getMRMOs()->getMaxMO());
RepresentationMatrix<MatrixType, VectorType> init2;

	RefSAFIndex = NULL;
	
	isSlave = 0;
}
*/




template <class MatrixType, class VectorType>
MRCIMatrix<MatrixType, VectorType>::~MRCIMatrix()
{
	if ( RefSAFIndex )
		delete[] RefSAFIndex;

	delete cache;

}


template <class MatrixType, class VectorType>
const NExternalsDiag *
	MRCIMatrix<MatrixType, VectorType>::getNExternalsDiag() const
{
	return mrccA;
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::_getDiagonal(MatrixType *p)
{
DiffConf<MOType>	diffConf;
TimeTicks	ticks;

	
	//	not very efficient
	ticks.start();
	for ( MRTreeIterator i = mrccA->firstInTree() ; 
			!mrccA->isLastInTree(i) ; mrccA->nextInTree(i) )
	{
		ConfigurationSAFNr<MOType> main = mrccA->getConfigurationSAFNr(i);
		diffConf.calcDiffConf(main, main);
		this->tablecase->calcLess3(diffConf);

		HMatElements<MatrixType, VectorType>	*repMats = 
			(HMatElements<MatrixType, VectorType> *)
			(*cache)[TableKey(*this->tablecase)];


	INT	safs = main.getSAFInc();
	MatrixType* pp = new MatrixType[safs*safs];
		repMats->getMatrix(pp, diffConf, *this->tablecase);

		for ( INT j=0 ; j<safs ; j++ )
			*p++ = pp[j*(safs+1)];
	delete[] pp;
	}
	ticks.stop();
	cout << ticks << endl;
}






template <class MatrixType, class VectorType>
LONG_LONG_INT MRCIMatrix<MatrixType, VectorType>::estimateNonZeroElements(
	INT maxBreak) const
{
LONG_LONG_INT	n = 0;

//	main loops over intern-ia / intern-ib
	for ( INT ia=0 ; ia<mrccA->getNumberOfElements() ; ia++ )
	{
		InternalConfsDiag	*mainContA = (*mrccA)[ia];

//	calculate upper triangular matrix only
		for ( INT ib=ia ; ib<mrccA->getNumberOfElements() ; ib++ )
		{
//	no interaction if difference in more than 2 MOs
			if ( abs(ia-ib)>2 )
				continue;

			InternalConfsDiag	*mainContB = (*mrccA)[ib];

			for ( INT ja=0 ; ja<mainContA->getNumberOfElements() ; ja++ )
			{
			TupelStructureDiag	*mainA = (*mainContA)[ja];		

//	calculate upper triangular matrix only
				for ( INT jb=(ia==ib) ? ja : 0 ; jb<mainContB->getNumberOfElements() ; jb++ )
				{
//					calculate cases with four external integrals seperately
//					if ( ia==2 && ib==2 && ja==jb )
//						continue;
						
					TupelStructureDiag	*mainB = (*mainContB)[jb];		

//		no interaction if more than double excitation between the two main-n
			INT	orderInt = 0;
					if ( (orderInt=Configuration<MOType>::calcExcitationOrderFast(*mainA, *mainB))>2 )
						continue;

				INT	maxExcitation = orderInt - (ib-ia) + ib;

					for ( INT ka=0 ; ka<mainA->getNumberOfElements() ; ka++ )
					{
					extMOsDiag	*MOA = (*mainA)[ka];
						for ( INT la=0 ; la<MOA->getNumberOfElements() ; la++ )
						{
						INT	jDiag = ib==ia && jb==ja;
							for ( INT kb=(jDiag) ? ka : 0 ; kb<mainB->getNumberOfElements() ; kb++ )
							{
							extMOsDiag	*MOB = (*mainB)[kb];
								if ( maxExcitation<=2 )
								{
									if ( jDiag )
										n += MOB->getNumberOfElements() *
											(MOA->getSAFInc()+1) * MOA->getSAFInc() / 2;
									else
										n += MOB->getNumberOfElements() *
											MOA->getSAFInc() * MOB->getSAFInc();
								}
								else
								if ( maxExcitation==3 )
										n += ((INT) sqrt((double)MOB->getNumberOfElements())) *
											MOA->getSAFInc() * MOB->getSAFInc();
								else
								if ( maxExcitation==4 )
										n += MOA->getSAFInc() * MOB->getSAFInc();
											
							if ( n>maxBreak )
								return n;
							}
						} // external MOs A
					} // tupel structure A
				} // internal-conf B
			} // internal-conf A
		} // internal-n B
	} // internal-n A
	
	return n;
}




template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::calcHamilton(
	NExternalsDiag	*mrccA, NExternalsDiag *mrccB, MultData data)
{
DiffConf<MOType>	diffConfExtAInternalB, diffConfInternal;
TimeTicks	ticks, ticksTotal;
TimeTicks	ticksSum;


INT	SAFsA = mrccA->getNumberOfTotalSpinAdaptedFunctions();

HMatElements<MatrixType, VectorType>	*hmat = NULL;

	switch ( data.calcMode ) {
	case MultData::_Multiplication:
		{
			ActionMultiplication<MatrixType, VectorType>	action(*data.x, *data.y);
			TwoDimIterator<
				MatrixType,
				VectorType, 
				ActionMultiplication<MatrixType, VectorType>,
				HMatElements<MatrixType, VectorType> >	iter(
					numberOfSlaves, &this->mrmos, this->Int4sharedMem, mrccA, mrccB, action, *hmat, 2, 1);
		}
		break;

	case MultData::_Density:
		{
		INT	nDens = 0;
		INT	n = data.x->getN();
			if ( data.y )
				n *= data.y->getN();
			for ( INT ix=0 ; ix<n ; ix++ )
			{
				memset(data.densMat[nDens++], 0, 
					mrccA->getMRMOs()->getMaxMO() * 
					mrccA->getMRMOs()->getMaxMO() * sizeof(VectorType));
			}

			ActionDensity<MatrixType, VectorType>	action(*data.x, *data.y, data.densMat);
			TwoDimIterator<
				MatrixType,
				VectorType, 
				ActionDensity<MatrixType, VectorType>,
				HMatElements<MatrixType, VectorType> >	iter(
					0*numberOfSlaves, &this->mrmos, this->Int4sharedMem, mrccA, mrccB, action, *hmat, 1, 0);
		}
		break;

	case MultData::_MatrixStorage:
		{
			ActionMatrixStorage<MatrixType, VectorType>	action(data.stored);
			TwoDimIterator<
				MatrixType,
				VectorType, 
				ActionMatrixStorage<MatrixType, VectorType>,
				HMatElements<MatrixType, VectorType> >	iter(
					0*numberOfSlaves, &this->mrmos, this->Int4sharedMem, mrccA, mrccB, action, *hmat, 2, 1);
		}
		break;

	case MultData::_PointerStorage:
		{
			memset(data.HMat, 0, SAFsA*mrccB->getNumberOfTotalSpinAdaptedFunctions()*sizeof(MatrixType));
			ActionPointerStorage<MatrixType, VectorType>	action(data.HMat, SAFsA);
			TwoDimIterator<
				MatrixType,
				VectorType, 
				ActionPointerStorage<MatrixType, VectorType>,
				HMatElements<MatrixType, VectorType> >	iter(
					0*numberOfSlaves, &this->mrmos, this->Int4sharedMem, mrccA, mrccB, action, *hmat, 2, 1);
		}
		break;
	}
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::_mult(DiskBuffer *xBuf, DiskBuffer *yBuf, INT start, INT end)
{
INT	freeMemory;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	freeMemory = 800000000;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
INT	Size = this->totalDim*sizeof(VectorType);

INT	perStep = freeMemory / (2*Size);
	if ( !perStep )
	{
		cout << "not enough memory, sorry" << endl;
		cerr << "not enough memory, sorry" << endl;
		exit(1);
	}
	if ( perStep>end - start + 1 )
		perStep = end - start + 1;


INT	steps = (end - start)/perStep + 1;


//	if ( getNumberOfSlaves() == 0 )
	{
	CIVectors<VectorType>	*x = new CIVectors<VectorType>(this->totalDim, perStep);
	CIVectors<VectorType>	*y = new CIVectors<VectorType>(this->totalDim, perStep);
//		cout << "performing ci-matrix generation and multiplication in " << steps;
//		cout << " step(s)" << endl;
//		cout << "with " << perStep << " vectors per step." << endl;
		
		for ( INT i=0 ; i<steps ; i++ )
		{	for ( INT j=0 ; j<perStep ; j++ )
				xBuf->get(start + i*perStep + j, x->getP(j));

			y->clear();



			calcHamilton(mrccA, mrccA, MultData(x, y));
			
	//		for ( INT ii=0 ; ii<this->totalDim ; ii++ )
	//			(*y)(ii) = rand();

			for ( INT j=0 ; j<perStep ; j++ )
				yBuf->put(start + i*perStep + j, y->getP(j));
//			cout << "x=" << *x << endl;
//			cout << "y=" << *y << endl;
//			cout << "step " << i+1 << " of " << steps << " finished." << endl;
		}
		delete x;
		delete y;
	}
/*
	else
	{
//		cout << "performing ci-matrix generation and multiplication with" << endl;
//		cout << getNumberOfSlaves() << " processes." << endl;
		for ( INT i=0 ; i<steps ; i++ )
		{	for ( INT j=0 ; j<perStep ; j++ )
				xBuf->get(start + i*perStep + j, getXVector()->getP(j));

			getYVector()->clear();

			calcHamilton();

			for ( INT j=0 ; j<perStep ; j++ )
				yBuf->put(start + i*perStep + j, getYVector()->getP(j));
//			cout << "y=" << *getYVector() << endl;
//			cout << "step " << i+1 << " of " << steps << " finished." << endl;
		}

	}
*/
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::getRefMatrix(MatrixType *pRefMatrix) const
{
NExternalsDiag	ref(mrccA->projectOnReferences());

	((MRCIMatrix *) this)->calcHamilton(&ref, &ref, MultData(pRefMatrix));
}

template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::getDensity(CIVectors<VectorType> & x, VectorType *densMat) const
{
	((MRCIMatrix *) this)->calcHamilton(mrccA, mrccA, MultData(&x, &x, &densMat));
}

template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::getReferenceDensity(CIVectors<VectorType> & x, VectorType *densMat) const
{
NExternalsDiag	ref(mrccA->projectOnReferences());

	((MRCIMatrix *) this)->calcHamilton(&ref, &ref, MultData(&x, &x, &densMat));
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::getTotalMatrix(MatrixStorage<MatrixType, VectorType> &m) const
{
	((MRCIMatrix *) this)->calcHamilton(mrccA, mrccA, MultData(&m));
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::getTotalMatrix(MatrixType *pCIMatrix) const
{
	((MRCIMatrix *) this)->calcHamilton(mrccA, mrccA, MultData(pCIMatrix));
}


template <class MatrixType, class VectorType>
void	MRCIMatrix<MatrixType, VectorType>::dens(
	const CIVectors<VectorType> * CIx, const CIVectors<VectorType> * CIy,
	VectorType **densMat) const
{
	((MRCIMatrix *) this)->calcHamilton(mrccA, mrccB, MultData(CIx, CIy, densMat));
}	


//#include "parallel.cch"


template class MRCIMatrix<float, float>;
template class MRCIMatrix<double, float>;
template class MRCIMatrix<float, double>;
template class MRCIMatrix<double, double>;


