//***********************************************************************
//
//	Name:			MRMPH0Matrix.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.09.1998
//
//
//
//
//
//***********************************************************************

#include "MRMPH0Matrix.h"

#include "../MRTree/Diag/NExternalsDiag.h"
#include "../MRTree/Diag/InternalConfsDiag.h"
#include "../MRTree/Diag/TupelStructureDiag.h"
#include "../MRTree/Diag/extMOsDiag.h"
#include "../MRTree/Container/ContainerIterator.h"

#include "../../Math/MatrixVector/BufferedVector.h"

#include "SelIterator.h"

#include "../Configuration/ConfigurationSet.h"

#include "../Cache/TableKey.h"
#include "../Cache/VarSizeReadOnlyCache.h"

#include "../MO/Iterators/MOIterator.h"
#include "../IO/TimeTicks.h"

#include "../MRTree/Diag/NExternalsDiag.Tmpl.h"
#include "MatrixAction/ActionE0.h"
#include "MatrixAction/ActionE1.h"

#include "../MO/Iterators/MOAccess.h"


#include <stdlib.h>
#include <iomanip>

using namespace std;

template <class MatrixType, class VectorType>
MRMPH0Matrix<MatrixType, VectorType>::MRMPH0Matrix(
	const CICalculation<MatrixType, VectorType> &ciCalc,
	const NExternalsDiag & Psi0,
	const NExternalsDiag & inhomogenity,
	const BufferedVector<VectorType> &Psi0EV,
	const BufferedVector<VectorType> &inhomogenityEV,
	const MRFockMatrix<VectorType> *mrFockMatrix,
	ProjectionMode projectionMode) :
		CICalculation<MatrixType, VectorType>(ciCalc), 
		Psi0(Psi0),
		inhomogenity(inhomogenity),
		Psi0EV(Psi0EV),
		inhomogenityEV(inhomogenityEV),
		alpha(NULL),
		projectionMode(projectionMode)
{
	init(mrFockMatrix);
}



template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::init(const MRFockMatrix<VectorType> *mrFockMatrix)
{
	irrep = Psi0.getTotalSymmetry();

	Hcache = new VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >
			(2048, (1<<20));
	H0cache = new VarSizeReadOnlyCache<TableKey, MRMPH0MatElements<MatrixType, VectorType> >
			(2048, (1<<20));

	buildInternalExcitations();
	createInteractionLists();

	E0 = 0;

MRMPH0MatElements<MatrixType, VectorType> init1(
	mrFockMatrix,
	HMatElements<MatrixType, VectorType>::int2Container,
	HMatElements<MatrixType, VectorType>::int4Container,
	HMatElements<MatrixType, VectorType>::int2Container->getCore(),
	E0);

	calcE0();
	calcE1();

MRMPH0MatElements<MatrixType, VectorType> init2(
	mrFockMatrix,
	HMatElements<MatrixType, VectorType>::int2Container,
	HMatElements<MatrixType, VectorType>::int4Container,
	HMatElements<MatrixType, VectorType>::int2Container->getCore(),
	E0);

	for ( INT i=0 ; i<=maxExc ; i++ )
		moAccess[i] = new MOAccess(i, &this->mrmos, MAXOPENSHELLS-2, Psi0);



	calcDim();


	if ( projectionMode == P_0Complement )
		calcAlpha();

/*	cout << "printing mat" << endl;
	printMat();
	cout << "OK" << endl;
*/
}



template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::calcE0()
{
TimeTicks	ticks;

INT	safs = Psi0.getNumberOfTotalSpinAdaptedFunctions();
VectorType	*Psi0EV1 = new VectorType[safs];

MatrixType	E0Tmp = 0;


	ticks.start();

	Psi0EV.get(Psi0EV1);

MRMPH0MatElements<MatrixType, VectorType>	*hmat = NULL;
ActionE0<MatrixType, VectorType>	action(Psi0EV1, E0Tmp);

	TwoDimIterator<
		MatrixType,
		VectorType, 
		ActionE0<MatrixType, VectorType>,
		MRMPH0MatElements<MatrixType, VectorType> >	iter(
			0, &this->mrmos, this->Int4sharedMem, &Psi0, &Psi0, action, *hmat, 1, 0);

	ticks.stop();
	cout << ticks << endl;

	delete Psi0EV1;

	E0 = E0Tmp;
	cout << "E0:" << E0 << endl;
	
}




template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::calcE1()
{
TimeTicks	ticks;

INT	safs = Psi0.getNumberOfTotalSpinAdaptedFunctions();
VectorType	*Psi0EV1 = new VectorType[safs];

	E1 = 0;

	ticks.start();

	Psi0EV.get(Psi0EV1);

HMatElements<MatrixType, VectorType>	*hmat = NULL;
ActionE1<MatrixType, VectorType>	action(Psi0EV1, E1);

	TwoDimIterator<
		MatrixType,
		VectorType, 
		ActionE1<MatrixType, VectorType>,
		HMatElements<MatrixType, VectorType> >	iter(
			0, &this->mrmos, this->Int4sharedMem, &Psi0, &Psi0, action, *hmat, 2, 0);

	ticks.stop();
	cout << ticks << endl;


	delete Psi0EV1;

// 	E1 -= E0[0];
	cout << "E0+E1:" << E1 << endl;
	
}





template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::calcDim()
{
	cout << "scanning MR-CI space... " << flush;


TimeTicks	ticks;

	ticks.start();
	this->totalDim = 0;
	totalConfs = 0;
	// use internal creators
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
		MOIterator	externalCreators(i, &this->mrmos, 1, 
			this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
			
			internal[i](iInt).setConfStart(totalConfs);
			internal[i](iInt).setSAFStart(this->totalDim);
			
		Configuration<MOType>	conf(internal[i](iInt));
		INT	open = conf.getNumberOfOpenShells();
			while ( !externalCreators.isEnd() )
			{
				totalConfs++;
				this->totalDim += Psi0.getNumberOfSpinAdaptedFunctions(
								open + externalCreators.getNumberOfOpenShells());
							
				externalCreators.next();
			}

			internal[i].next(iInt);
		}
	}
	ticks.stop();
	cout << ticks << endl;
	cout << "totalDim=" << this->totalDim << endl;
	cout << "totalConfs=" << totalConfs << endl;


/*	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "::::::::::::::::: i=" << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			cout << ((Configuration<MOType>) internal[i](iInt)) << "  " << 
				internal[i](iInt).getSAFStart() << " " << 
				internal[i](iInt).calcIrRep(this->mrmos) << endl;
			internal[i].next(iInt);
		}
	}
*/
}


template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::buildInternalExcitations()
{
	cout << "creating internal space... " << flush;

TimeTicks	ticks;
	ticks.start();

	

ConfigurationSet	references(Psi0.getRefConfSet());

/*
	for ( INT ia=0 ; ia<Psi0.getNumberOfElements() ; ia++ )
	{
	const InternalConfsDiag	*mainContA = Psi0[ia];

		for ( INT ja=0 ; ja<mainContA->getNumberOfElements() ; ja++ )
		{
		const TupelStructureDiag	*mainA = (*mainContA)[ja];
			if ( mainA->isReference() )
			{
			Configuration<MOType>	conf(*mainA);
				references.add(conf);
			}
		}
	}
*/
	

	// use internal annihilators
Pix	iRef = references.first();
	while ( iRef )
	{
		for ( INT i=0 ; i<=maxExc ; i++ )
		{
		MOIterator	internalAnnihilators(i, &this->mrmos, 0);
			while ( !internalAnnihilators.isEnd() )
			{
			ConfigurationStart<MOType>	conf(references(iRef));
				conf -= internalAnnihilators;
				if ( conf.getNumberOfElectrons() )
					internal[i].add(conf);
				internalAnnihilators.next();
			}
		}
		references.next(iRef);
	}

/*	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "::::::::::::::::: i=" << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			cout << ((Configuration<MOType>) internal[i](iInt)) << "  " << 
				internal[i](iInt).calcIrRep(this->mrmos) << endl;
			internal[i].next(iInt);
		}
	}
	cout << "-------------------------------------" << endl;
*/

	// use internal creators
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			for ( INT j=1 ; j<=i ; j++ )
			{
			MOIterator	*internalCreators;

				if ( i-j==0 )
					internalCreators = new MOIterator(j, &this->mrmos, 0,
						this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
				else
					internalCreators = new MOIterator(j, &this->mrmos, 0);

				while ( !internalCreators->isEnd() )
				{
				ConfigurationStart<MOType>	conf(internal[i](iInt));
//					cout << ((Configuration<MOType>) conf) << " " << *internalCreators << " ";
					conf += *internalCreators;
//					cout << ((Configuration<MOType>) conf) << endl;

					if ( conf.getNumberOfElectrons() )
						internal[i-j].add(conf);
					internalCreators->next();
				}

				delete internalCreators;
			}
			internal[i].next(iInt);
		}
	}


/*	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "::::::::::::::::: i=" << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			cout << ((Configuration<MOType>) internal[i](iInt)) << "  " << 
				internal[i](iInt).getSAFStart() << " " << 
				internal[i](iInt).calcIrRep(this->mrmos) << endl;
			internal[i].next(iInt);
		}
	}
*/

	ticks.stop();
	cout << ticks << endl;
}



template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::createInteractionLists()
{
	cout << "creating internal interaction lists... " << flush;

TimeTicks	ticks;
	ticks.start();


	interactionList = new Pix *** [maxExc+1];
	interactionSymList = new IrRep *** [maxExc+1];
	nInteractions = new INT ** [maxExc+1];
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		interactionList[i] = new Pix ** [internal[i].length()];
		interactionSymList[i] = new IrRep ** [internal[i].length()];
		nInteractions[i] = new INT * [internal[i].length()];
	Pix	iInt = internal[i].first();
	INT	iInti = 0;
		while ( iInt )
		{
			interactionList[i][iInti] = new Pix * [maxExc+1];
			interactionSymList[i][iInti] = new IrRep * [maxExc+1];
			nInteractions[i][iInti] = new INT [maxExc+1];
		Configuration<MOType>	iConf(internal[i](iInt));

//			cout << iConf << endl;
			for ( INT j=i ; j<=maxExc ; j++ )
			{
				interactionList[i][iInti][j] = NULL;
				interactionSymList[i][iInti][j] = 0;
				nInteractions[i][iInti][j] = 0;
			
				// no interaction if more than 1 excitation
				if ( abs(i-j)>1 )
					continue;


				// scan number of interactions
			INT* exc = new INT[internal[j].length()];
				memset(exc, 0, internal[j].length()*sizeof(INT));
			Pix	jInt = (i==j ? iInt : internal[j].first());
			INT	jInti = 0;
			INT	n = 0;
				if ( j>=i )
					while ( jInt )
					{
						if ( Configuration<MOType>::calcExcitationOrderFast(
							iConf, internal[j](jInt), 1)<=1 )
						{
							exc[jInti] = 1;
							n++;
						}
						internal[j].next(jInt);
						jInti++;
					}
				else
					while ( jInt )
					{
						if ( Configuration<MOType>::calcExcitationOrderFast(
							internal[j](jInt), iConf, 1)<=1 )
						{
							exc[jInti] = 1;
							n++;
						}
						internal[j].next(jInt);
						jInti++;
					}



				// store interactions
				interactionList[i][iInti][j] = new Pix [n];
				interactionSymList[i][iInti][j] = new IrRep [n];
				nInteractions[i][iInti][j] = n;
				jInt = (i==j ? iInt : internal[j].first());
				jInti = 0;
			Pix	*p = interactionList[i][iInti][j];
			IrRep	*pSym = interactionSymList[i][iInti][j];
				while ( jInt )
				{
					if ( exc[jInti] )
					{
						*p++ = jInt;
						*pSym++ = internal[j](jInt).calcIrRep(this->mrmos);
					}

					internal[j].next(jInt);
					jInti++;
				}
				delete[] exc;

			}
			internal[i].next(iInt);
			iInti++;
		}
	}

	ticks.stop();
	cout << ticks << endl;
}



template <class MatrixType, class VectorType>
MRMPH0Matrix<MatrixType, VectorType>::~MRMPH0Matrix()
{
	delete Hcache;
	delete H0cache;


	if ( alpha ) 
		delete alpha;
		
		
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		delete moAccess[i];

	Pix	iInt = internal[i].first();
	INT	iInti = 0;
		while ( iInt )
		{

			for ( INT j=i ; j<=maxExc ; j++ )
			{
				if ( interactionList[i][iInti][j] )
					delete interactionList[i][iInti][j];
				if ( interactionSymList[i][iInti][j] )
					delete interactionSymList[i][iInti][j];
			}

			if ( interactionList[i][iInti] )
				delete interactionList[i][iInti];
			if ( interactionSymList[i][iInti] )
				delete interactionSymList[i][iInti];
			if ( nInteractions[i][iInti] )
				delete nInteractions[i][iInti];

			internal[i].next(iInt);
			iInti++;
		}
		if ( interactionList[i] )
			delete interactionList[i];
		if ( interactionSymList[i] )
			delete interactionSymList[i];
		if ( nInteractions[i] )
			delete nInteractions[i];
	}
	if ( interactionList )
		delete interactionList;
	if ( interactionSymList )
		delete interactionSymList;
	if ( nInteractions )
		delete nInteractions;
}


template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::clearPsi0(BufferedVector<VectorType> &x) const
{
TimeTicks	ticks;


	cout << "zeroing Psi0" << endl;


	ticks.start();

SelIterator	sel(Psi0, internal, moAccess, maxExc);

	while ( !sel.isEnd() )
	{
	INT	n = sel.getSAFStart();
		for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
			x[n+i] = 0;
		sel.next();
	}

	ticks.stop();
	cout << ticks << endl;
}


template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::printMat() const
{
DiffConf<MOType>	diffConf;

MatrixType	*mat = new MatrixType[this->totalDim*this->totalDim];
	memset(mat, 0, this->totalDim*this->totalDim*sizeof(MatrixType));

	// use internal creators
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
	Pix	iInt = internal[i].first();
		while ( iInt )
		{

		INT	iSafStart = internal[i](iInt).getSAFStart();

		MOIterator	iexternalCreators(i, &this->mrmos, 1, 
			this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));

			while ( !iexternalCreators.isEnd() )
			{
			Configuration<MOType>	iconf(internal[i](iInt));
				iconf += iexternalCreators;

			INT	iSafs = Psi0.getNumberOfSpinAdaptedFunctions(
				iconf.getNumberOfOpenShells());


				for ( INT j=0 ; j<=maxExc ; j++ )
				{
				Pix	jInt = internal[j].first();
					while ( jInt )
					{

					INT	jSafStart = internal[j](jInt).getSAFStart();

					MOIterator	jexternalCreators(j, &this->mrmos, 1, 
						this->mrmos.getProd(internal[j](jInt).calcIrRep(this->mrmos), irrep));

						while ( !jexternalCreators.isEnd() )
						{
						Configuration<MOType>	jconf(internal[j](jInt));
							jconf += jexternalCreators;

						INT	jSafs = Psi0.getNumberOfSpinAdaptedFunctions(
							jconf.getNumberOfOpenShells());



							diffConf.calcDiffConf(iconf, jconf);
							this->tablecase->calcLess3(diffConf);

							const MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
								(*H0cache)[TableKey(*this->tablecase)];
//								const HMatElements<MatrixType, VectorType>	*repMats = 
//									(*Hcache)[TableKey(*this->tablecase)];



						MatrixType* pp = new MatrixType[iSafs*jSafs];

						MatrixType	*pMat = pp;

//								repMats->getMatrix(pp, diffConf, *this->tablecase);
							repMats->getMatrix(pp, diffConf, *this->tablecase);

							for ( INT i=0 ; i<iSafs ; i++ )
								for ( INT j=0 ; j<jSafs ; j++ )
									mat[(iSafStart+i)*this->totalDim + jSafStart + j] = *pMat++;


							jSafStart += jSafs;
							jexternalCreators.next();
							delete[] pp;
						}
						internal[j].next(jInt);
					}
				}
				iSafStart += iSafs;
				iexternalCreators.next();
			}
			internal[i].next(iInt);
		}
	}
	
	for ( INT i=0 ; i<this->totalDim ; i++ )
	{
		for ( INT j=0 ; j<this->totalDim ; j++ )
		{
			cout << setw(14) << mat[i*this->totalDim + j] << "\t";
		}
		cout << endl;
	}

	delete mat;
}



template <class MatrixType, class VectorType>
INT	MRMPH0Matrix<MatrixType, VectorType>::getConfNr(const Configuration<MOType> &conf) const
{
ConfigurationStart<MOType>	intPart;
Configuration<MOType>	extPart;
	conf.split(&this->mrmos, intPart, extPart);
	
INT	nExt = extPart.getNumberOfElectrons();
Pix	pix = internal[nExt].seek(intPart);
	if ( !pix )
	{
		cout << "wrong Psi0 wave function" << endl;
		cout << ((Configuration<MOType>) conf) << endl;
		exit(1);
	}

	INT	ConfStart = internal[nExt](pix).getConfStart();

	switch ( nExt ) {
	case 0:
		return ConfStart;

	case 1:
		return ConfStart + moAccess[nExt]->getConfOffset(extPart.getOpenShell(0));

	case 2:
		if ( extPart.getNumberOfOpenShells()==2 )
			return ConfStart + moAccess[nExt]->getConfOffset(
				extPart.getOpenShell(0), extPart.getOpenShell(1));
		else
			return ConfStart + moAccess[nExt]->getConfOffset(
				extPart.getClosedShell(0), extPart.getClosedShell(0));
	}
	return 0;
}



template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::calcResidual(BufferedVector<VectorType> &x) const
// calculate resiual vector
//          x  = b - x
{
DiffConf<MOType>	diffConf;
TimeTicks	ticks;

	cout << "performing x = b - x" << endl;
	

	ticks.start();



//	switch ( inhomogenity ) {
//	case I_Psi0:
		{
		INT	safs = inhomogenity.getNumberOfTotalSpinAdaptedFunctions();
		VectorType	*inhomogenityEV1 = new VectorType[safs];

			inhomogenityEV.get(inhomogenityEV1);
		Configuration<MOType>	*from = (Configuration<MOType>	*) &diffConf.getFrom();
			for ( INT ia=0 ; ia<inhomogenity.getNumberOfElements() ; ia++ )
			{
			const InternalConfsDiag	*mainContA = inhomogenity[ia];

				for ( INT ja=0 ; ja<mainContA->getNumberOfElements() ; ja++ )
				{
				const TupelStructureDiag	*mainA = (*mainContA)[ja];		


					for ( INT i=0 ; i<=maxExc ; i++ )
					{
					Pix	iInt = internal[i].first();
						while ( iInt )
						{
						INT	internalOrder = Configuration<MOType>::calcExcitationOrder(
								*mainA, internal[i](iInt));

							if ( internalOrder>2 )
							{
								internal[i].next(iInt);
								continue;
							}


							for ( INT ka=0 ; ka<mainA->getNumberOfElements() ; ka++ )
							{
							const extMOsDiag	*MOs = (*mainA)[ka];
								if ( !MOs )
									continue;

							INT	iSafStart = MOs->getSAFStart();
							INT	iSafs = MOs->getSAFInc();

		//	cout << "Hallo1 " << MOs << " " << MOs->getNumberOfTotalMOs() << "... " << iSafs << endl;
							const MOType	*mo = NULL;
								if ( MOs->getNumberOfTotalMOs() && MOs->getNumberOfElements() )
									mo = (*((extMOsDiag *) MOs))[0];
		//	cout << "Hallo2" << endl;

								for ( INT la=0 ; la<MOs->getNumberOfElements() ; la++ )
								{
								Configuration<MOType> ref = (*mainA);
		//								cout << "A" << ref << endl;

								MOType	extA[2] = { 0, 0 };

									for ( INT k=0 ; k<MOs->getNumberOfOpenMOs() ; k++ )
									{
										extA[k] = *mo;
										ref.insertOpenMO(*mo++);
									}

									if ( MOs->getNumberOfClosedMOs() )
									{
										extA[0] = extA[1] = *mo;
										ref.insertClosedMO(*mo++);
									}



		//								cout << "B" << iSafStart << endl; //ref << endl;

		//								cout << ((Configuration<MOType> &) *mainA) << ", SAFStart= " << iSafStart << endl;

							INT	xOffset = internal[i](iInt).getSAFStart();
							INT	jSafStart = internal[i](iInt).getSAFStart();
							INT	internalOpenB = internal[i](iInt).getNumberOfOpenShells();
							INT	CmpStatusOld = -1;
							INT	QCmpStatusOld = 0;
							MatrixType* mat = new MatrixType[iSafs*inhomogenity.getNumberOfSpinAdaptedFunctions(
											internalOpenB+2)];
							INT	extPosBOpen[2] = { -1, -1 };
							INT	extPosBClosed[2] = { -1, -1 };

								MOIterator	externalCreators(i, &this->mrmos, 1, 
									this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
									while ( !externalCreators.isEnd() )
									{

								MOType	extB[2] = { this->mrmos.getMaxMO()+1, this->mrmos.getMaxMO()+2 };

										if ( i>0 )
										{
											extB[0] = externalCreators.getMO(0);
											if ( i>1 )
												extB[1] = externalCreators.getMO(1);
										}

									INT	order = (i<ia ? i : ia);
										if ( extA[0]==extB[0] )
										{
											order--;
											order -= extA[1]==extB[1];
										}
										else
										{
											if ( extA[1]==extB[1] )
												order--;
											else
											{
												order -= extA[1]==extB[0];
												order -= extA[0]==extB[1];
											}
										}

									INT	open = i - ((extB[0]==extB[1]) << 1);

										if ( internalOrder + order > 2 )
										{
											xOffset += inhomogenity.getNumberOfSpinAdaptedFunctions(
												internalOpenB + open);
											externalCreators.next();
											continue;
										}


									INT	CmpStatus = 
											((extA[0]==extB[0]) << 0) |
											((extA[1]==extB[1]) << 1) |
											((extA[0]==extB[1]) << 2) |
											((extA[1]==extB[0]) << 3) |
											((extB[0]==extB[1]) << 4);

									INT	QCmpStatus = 
											((extA[0] <extB[0]) << 0) |
											((extA[1] <extB[1]) << 1) |
											((extA[0] <extB[1]) << 2) |
											((extA[1] <extB[0]) << 3);


										if ( externalCreators.changedCase() || 
											CmpStatus!=CmpStatusOld || QCmpStatus!=QCmpStatusOld )
										{
										Configuration<MOType>	conf(internal[i](iInt));
											conf += externalCreators;

		//										cout << "conf=" << conf << endl;
											diffConf.calcDiffConf(conf, ref);


											extPosBOpen[0] = extPosBOpen[1] = extPosBClosed[0] = extPosBClosed[1] = -1;

											for ( INT ii=0 ; ii<from->getNumberOfOpenShells() ; ii++ )
											{
												if ( extB[0] == from->getOpenShell(ii) )
													extPosBOpen[0] = ii;
												if ( extB[1] == from->getOpenShell(ii) )
													extPosBOpen[1] = ii;
											}

											for ( INT ii=0 ; ii<from->getNumberOfClosedShells() ; ii++ )
											{
												if ( extB[0] == from->getClosedShell(ii) )
													extPosBClosed[0] = ii;
												if ( extB[1] == from->getClosedShell(ii) )
													extPosBClosed[1] = ii;
											}

										}
										CmpStatusOld = CmpStatus;
										QCmpStatusOld = QCmpStatus;



		//									cout << diffConf << endl;
										//	update MOs in DiffConf
										if ( extPosBOpen[0]>=0 )
											from->setOpenShell(extPosBOpen[0], extB[0]);
										if ( extPosBOpen[1]>=0 )
											from->setOpenShell(extPosBOpen[1], extB[1]);

										if ( extPosBClosed[0]>=0 )
											from->setClosedShell(extPosBClosed[0], extB[0]);
										if ( extPosBClosed[1]>=0 )
											from->setClosedShell(extPosBClosed[1], extB[1]);

		//									cout << diffConf << endl;
										this->tablecase->calcLess3(diffConf);
		//									cout << *this->tablecase << endl;

										const HMatElements<MatrixType, VectorType>	*repMats = 
											(*Hcache)[TableKey(*this->tablecase)];


										jSafStart += repMats->getNumberOfRows();


										repMats->getMatrix(mat, diffConf, *this->tablecase);


									MatrixType	*pMat = mat;
										for ( INT j=0 ; j<repMats->getNumberOfRows() ; j++ )
										{
										VectorType	h = 0;
											for ( INT i=0 ; i<iSafs ; i++ )
											{
		//											cout << "i=" << i << ", j=" << j << ": " << 
		//												*pMat << " " << inhomogenityEV1[iSafStart+i] << " " <<
		//												*pMat * inhomogenityEV1[iSafStart+i] << endl;
												h += *pMat++ * inhomogenityEV1[iSafStart+i];
											}
		//										cout << "h=" << h << endl;
		//										cout << "offset=" << p-v << endl;
											x[xOffset++] += h;
										}



										externalCreators.next();
									}
									iSafStart += iSafs;
									delete[] mat;
								}
							}
							internal[i].next(iInt);
						}
					}
				}
			}
			delete inhomogenityEV1;
		}
/*		break;
		
	case I_Reference:
		{
		ConfigurationSet	references(Psi0.getRefConfSet());
		INT	nMax = Psi0.getNumberOfSpinAdaptedFunctions(MAXSGATABLEOPENSHELLS);
		MatrixType	mat[nMax*nMax];

		INT	safs = Psi0.getNumberOfRefConfSpinAdaptedFunctions();
		Vector<VectorType>	ev(safs);
		
			{
			INT	m1 = 0;
			INT	m2 = 0;
				for ( INT ia=0 ; ia<Psi0.getNumberOfElements() ; ia++ )
				{
				const InternalConfsDiag	*mainContA = Psi0[ia];

					for ( INT ja=0 ; ja<mainContA->getNumberOfElements() ; ja++ )
					{
					const TupelStructureDiag	*mainA = (*mainContA)[ja];
					INT	n = Psi0.getNumberOfSpinAdaptedFunctions(
							mainA->getNumberOfOpenShells());
							
						if ( mainA->isReference() )
						{
							for ( INT k=0 ; k<n ; k++ )
								ev[m1++] = Psi0EV[m2++];
						}
						else
							m2 += n;
					}
				}
			}

		
			ev.normalize2();
			
		
		INT	iSafStart = 0;

		Pix	iRef = references.first();
			while ( iRef )
			{
				ConfigurationStart<MOType>	ref(references(iRef));

				for ( INT i=0 ; i<=maxExc ; i++ )
				{
				Pix	iInt = internal[i].first();
					while ( iInt )
					{
					INT	internalOrder = Configuration<MOType>::calcExcitationOrder(
							ref, internal[i](iInt));

						if ( internalOrder>2 )
						{
							internal[i].next(iInt);
							continue;
						}
					INT	xOffset = internal[i](iInt).getSAFStart();
					MOIterator	externalCreators(i, &this->mrmos, 1, 
						this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
						while ( !externalCreators.isEnd() )
						{
						Configuration<MOType>	conf(internal[i](iInt));
							conf += externalCreators;

							diffConf.calcDiffConf(conf, ref);
							this->tablecase->calcLess3(diffConf);

							const HMatElements<MatrixType, VectorType>	*repMats = 
								(*Hcache)[TableKey(*this->tablecase)];

							repMats->getMatrix(mat, diffConf, *this->tablecase);

						MatrixType	*pMat = mat;
							for ( INT j=0 ; j<repMats->getNumberOfRows() ; j++ )
							{
							VectorType	h = 0;
								for ( INT i=0 ; i<repMats->getNumberOfColumns() ; i++ )
								{
									h += *pMat++ * ev[iSafStart+i];
								}
								x[xOffset++] += h;
							}
							externalCreators.next();
						}
						internal[i].next(iInt);
					}
				}
				iSafStart += Psi0.getNumberOfSpinAdaptedFunctions(
					ref.getNumberOfOpenShells());
				references.next(iRef);
			}
		}
		break;
	}
*/
	ticks.stop();
	cout << ticks << endl;


//	cout << "b-x=" << endl;
//	for ( INT i=0 ; i<this->totalDim ; i++ )
//		cout << v[i] << endl;	
	
	
	switch ( projectionMode ) {
	case P_no0:
		clearPsi0(x);
		break;
	case P_0Complement:
		//-------------------------------
		//	project out eigenvector
		{
		INT	safs = Psi0.getNumberOfTotalSpinAdaptedFunctions();
		VectorType	*ev = new VectorType[safs];
			Psi0EV.get(ev);

		double	evProj = 0;
			{
			INT m = 0;
			SelIterator	sel(Psi0, internal, moAccess, maxExc);
				while ( !sel.isEnd() )
				{
				INT	n = sel.getSAFStart();
					for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
						evProj += x[n++]*ev[m++];
					sel.next();
				}
			}

			cout << "evProj=" << evProj << endl;

			{
			INT m = 0;
			SelIterator	sel(Psi0, internal, moAccess, maxExc);
				while ( !sel.isEnd() )
				{
				INT	n = sel.getSAFStart();
					for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
						x[n++] -= evProj*ev[m++];
					sel.next();
				}
			}
		}
		//-------------------------------
		break;
	}
}



template <class MatrixType, class VectorType>
MatrixType	MRMPH0Matrix<MatrixType, VectorType>::calcMP3(const BufferedVector<VectorType> &x) const
{
DiffConf<MOType>	diffConf;
TimeTicks	ticks;

	cout << "performing E3 = <Psi|H1-E1|Psi>" << endl;
	
VectorType	*v = new VectorType[this->totalDim];
MatrixType	E3 = 0;
Configuration<MOType>	*from = (Configuration<MOType>	*) &diffConf.getFrom();

	ticks.start();

	x.get(v);


	for ( INT ia=0 ; ia<=maxExc ; ia++ )
	{
		for ( INT ib=ia ; ib<=maxExc ; ib++ )
		{
		Pix	iaInt = internal[ia].first();
			while ( iaInt )
			{
			Pix	ibInt = (ib==ia ? iaInt : internal[ib].first());
				while ( ibInt )
				{
				INT	internalOrder = Configuration<MOType>::calcExcitationOrder(
						internal[ia](iaInt), internal[ib](ibInt));

					if ( internalOrder>2 )
					{
						internal[ib].next(ibInt);
						continue;
					}


				VectorType	*pia = v + internal[ia](iaInt).getSAFStart();
				MOIterator	aExternalCreators(ia, &this->mrmos, 1, 
					this->mrmos.getProd(internal[ia](iaInt).calcIrRep(this->mrmos), irrep));
				INT	aExti = 0;
					while ( !aExternalCreators.isEnd() )
					{
					Configuration<MOType>	ref(internal[ia](iaInt));
						ref += aExternalCreators;

				MOType	extA[2] = { this->mrmos.getMaxMO()+1, this->mrmos.getMaxMO()+2 };

						if ( ia>0 )
						{
							extA[0] = aExternalCreators.getMO(0);
							if ( ia>1 )
								extA[1] = aExternalCreators.getMO(1);
						}

					VectorType	*pib = v + internal[ib](ibInt).getSAFStart();
					INT	internalOpenB = internal[ib](ibInt).getNumberOfOpenShells();

					INT	CmpStatusOld = -1;
					INT	QCmpStatusOld = 0;
					INT	aSAFInc = Psi0.getNumberOfSpinAdaptedFunctions(
										ref.getNumberOfOpenShells());
					MatrixType* mat = new MatrixType[aSAFInc*Psi0.getNumberOfSpinAdaptedFunctions(
										internalOpenB+2)];
					MatrixType* matH0 = new MatrixType[aSAFInc*Psi0.getNumberOfSpinAdaptedFunctions(
										internalOpenB+2)];
					INT	extPosBOpen[2] = { -1, -1 };
					INT	extPosBClosed[2] = { -1, -1 };



					MOIterator	bExternalCreators(ib, &this->mrmos, 1, 
						this->mrmos.getProd(internal[ib](ibInt).calcIrRep(this->mrmos), irrep));
						if ( ia==ib && iaInt==ibInt )
						{
						INT	bExti = 0;
							while ( !bExternalCreators.isEnd() && bExti<aExti )
							{
							INT	open = ib;
								if ( ib==2 && 
									bExternalCreators.getMO(0)==bExternalCreators.getMO(1) )
									open -= 2;
								pib += Psi0.getNumberOfSpinAdaptedFunctions(
									internalOpenB + open);
								bExternalCreators.next();
								bExti++;
							}
						}
						while ( !bExternalCreators.isEnd() )
						{

					MOType	extB[2] = { this->mrmos.getMaxMO()+3, this->mrmos.getMaxMO()+4 };

							if ( ib>0 )
							{
								extB[0] = bExternalCreators.getMO(0);
								if ( ib>1 )
									extB[1] = bExternalCreators.getMO(1);
							}



//------------------------------------------------------------------------------------

						INT	order = (ib<ia ? ib : ia);
							if ( extA[0]==extB[0] )
							{
								order--;
								order -= extA[1]==extB[1];
							}
							else
							{
								if ( extA[1]==extB[1] )
									order--;
								else
								{
									order -= extA[1]==extB[0];
									order -= extA[0]==extB[1];
								}
							}

						INT	open = ib - ((extB[0]==extB[1]) << 1);

							if ( internalOrder + order > 2 )
							{
								pib += Psi0.getNumberOfSpinAdaptedFunctions(
									internalOpenB + open);
								bExternalCreators.next();
								continue;
							}


						INT	CmpStatus = 
								((extA[0]==extB[0]) << 0) |
								((extA[1]==extB[1]) << 1) |
								((extA[0]==extB[1]) << 2) |
								((extA[1]==extB[0]) << 3) |
								((extB[0]==extB[1]) << 4);

						INT	QCmpStatus = 
								((extA[0] <extB[0]) << 0) |
								((extA[1] <extB[1]) << 1) |
								((extA[0] <extB[1]) << 2) |
								((extA[1] <extB[0]) << 3);


							if ( bExternalCreators.changedCase() || 
								CmpStatus!=CmpStatusOld || QCmpStatus!=QCmpStatusOld )
							{
							Configuration<MOType>	conf(internal[ib](ibInt));
								conf += bExternalCreators;

								diffConf.calcDiffConf(conf, ref);

								extPosBOpen[0] = extPosBOpen[1] = extPosBClosed[0] = extPosBClosed[1] = -1;

								for ( INT ii=0 ; ii<from->getNumberOfOpenShells() ; ii++ )
								{
									if ( extB[0] == from->getOpenShell(ii) )
										extPosBOpen[0] = ii;
									if ( extB[1] == from->getOpenShell(ii) )
										extPosBOpen[1] = ii;
								}

								for ( INT ii=0 ; ii<from->getNumberOfClosedShells() ; ii++ )
								{
									if ( extB[0] == from->getClosedShell(ii) )
										extPosBClosed[0] = ii;
									if ( extB[1] == from->getClosedShell(ii) )
										extPosBClosed[1] = ii;
								}

							}
							CmpStatusOld = CmpStatus;
							QCmpStatusOld = QCmpStatus;


							//	update MOs in DiffConf
							if ( extPosBOpen[0]>=0 )
								from->setOpenShell(extPosBOpen[0], extB[0]);
							if ( extPosBOpen[1]>=0 )
								from->setOpenShell(extPosBOpen[1], extB[1]);

							if ( extPosBClosed[0]>=0 )
								from->setClosedShell(extPosBClosed[0], extB[0]);
							if ( extPosBClosed[1]>=0 )
								from->setClosedShell(extPosBClosed[1], extB[1]);

							this->tablecase->calcLess3(diffConf);

						const HMatElements<MatrixType, VectorType>	*repMats = 
								(*Hcache)[TableKey(*this->tablecase)];

							repMats->getMatrix(mat, diffConf, *this->tablecase);


						const MRMPH0MatElements<MatrixType, VectorType>	*repMatsH0 = NULL;

							if ( this->tablecase->getP()==3 || this->tablecase->getP()==5 )
							{
								repMatsH0 = 
									(*H0cache)[TableKey(*this->tablecase)];
                                                                  repMatsH0->getMatrix(matH0, diffConf, *this->tablecase);
                                                                  for ( INT k=0 ; k<repMats->getNumberOfRows()*
                                                                                  repMats->getNumberOfColumns() ; k++ )
                                                                          mat[k] -= matH0[k];
                                                          }



                                                  MatrixType	*pMat = mat;

  /*								for ( INT jb=0 ; jb<repMats->getNumberOfRows() ; jb++ )
                                                          {
                                                                  for ( INT ja=0 ; ja<repMats->getNumberOfColumns() ; ja++ )
                                                                          cout << setw(16) << *pMat++ << "\t";
                                                                  cout << endl;
                                                          }
                                                          pMat = mat;
  */								
                                                  INT	diag = (pia==pib);
                                                          if ( diag )
                                                                  for ( INT jb=0 ; jb<repMats->getNumberOfRows() ; jb++ )
                                                                  {
                                                                          pMat += jb;
  //										cout << "jb=" << jb << " " << *pMat << " " << pib[jb] << " " << pia[jb] << endl;
                                                                          E3 += pib[jb] * (*pMat++ - E1)* pia[jb];
                                                                          for ( INT ja=jb+1 ; ja<repMats->getNumberOfColumns() ; ja++ )
                                                                                  E3 += 2 * pib[jb] * *pMat++ * pia[ja];
                                                                  }
                                                          else
                                                                  for ( INT jb=0 ; jb<repMats->getNumberOfRows() ; jb++ )
                                                                          for ( INT ja=0 ; ja<repMats->getNumberOfColumns() ; ja++ )
                                                                          {
  //											cout << "jb=" << jb << ", ja=" << ja << " " << *pMat << " " << pib[jb] << " " << pia[ja] << endl;
                                                                                  E3 += 2 * pib[jb] * *pMat++ * pia[ja];
  //											cout << "E3=" << E3 << endl;
                                                                          }

                                                          pib += repMats->getNumberOfRows();
  //------------------------------------------------------------------------------------
                                                          bExternalCreators.next();
                                                  }
                                                  pia += aSAFInc;
                                                  aExternalCreators.next();
                                                  aExti++;
                                                  delete[] matH0;
                                                  delete[] mat;
                                          }
                                          internal[ib].next(ibInt);
                                  }
                                  internal[ia].next(iaInt);
                          }
                  }
          }
          ticks.stop();
          cout << ticks << endl;
          
          delete v;
          
          return E3;
  }





  template <class MatrixType, class VectorType>
  void	MRMPH0Matrix<MatrixType, VectorType>::multInvDiag(BufferedVector<VectorType> &x) const
  // perform multiplication with invers of diagonal
  //                  1
  //          x  = -------  x 
  //           i      A      i
  //                   ii
  {
  TimeTicks	ticks;

  DiffConf<MOType>	diffConf;


          cout << "performing x = 1/D * x" << endl;

  INT	nn = 0;
          ticks.start();

  INT	xOffset = 0;
          // use internal creators
          for ( INT i=0 ; i<=maxExc ; i++ )
          {
          Pix	iInt = internal[i].first();
                  while ( iInt )
                  {
                  MOIterator	externalCreators(i, &this->mrmos, 1, 
                          this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
                          while ( !externalCreators.isEnd() )
                          {
                          Configuration<MOType>	conf(internal[i](iInt));
                                  conf += externalCreators;
                          INT	safs = Psi0.getNumberOfSpinAdaptedFunctions(
                                  conf.getNumberOfOpenShells());
                                  nn += safs;
  //					cout << "nn=" << nn << endl;


  //					cout << conf << ": ";

                          MatrixType	diag = MRMPH0MatElements<MatrixType, VectorType>::getDiag(conf);

                                  diag -= E0;

                                  for ( INT j=0 ; j<safs ; j++ )
                                          x[xOffset++] /= diag;

  //					cout << "diag=" << diag << ", E0=" << E0[r] << ", diff=" << diag-E0[r] << endl;



  /*				diffConf.calcDiffConf(conf, conf);
                          this->tablecase->calcLess3(diffConf);
		MatrixType	pp[safs*safs];

			const HMatElements<MatrixType, VectorType>	*repMats = 
				(*Hcache)[TableKey(*this->tablecase)];
			repMats->getMatrix(pp, diffConf, *this->tablecase);

				for ( INT j=0 ; j<safs ; j++ )
				{
//						cout << pp[j*(safs+1)] << " " << E0[r] << " " << pp[j*(safs+1)]+E0[r] << " " << *p/pp[j*(safs+1)] << endl;
					x[XOffset++] /= pp[j*(safs+1)];
				}

*/
				externalCreators.next();
			}
			internal[i].next(iInt);
		}
	}
	ticks.stop();
	cout << ticks << endl;
	
	switch ( projectionMode ) {
	case P_no0:
		clearPsi0(x);
		break;

	case P_0Complement:
		break;
	}
	
}


template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::multFock(
		const BufferedVector<VectorType> &x, 
		BufferedVector<VectorType> &y) const
// perform multiplication on n vectors
// y = A*x
{
TimeTicks	ticks;
DiffConf<MOType>	diffConf;


VectorType	*vx = new VectorType[this->totalDim];
VectorType	*vy = new VectorType[this->totalDim];
//	cout << "totalDim=" << this->totalDim << ", " << x->getObjectSize() << endl;
	cout.setf(ios::scientific);

	cout << "performing y = A*x" << endl;

	ticks.start();

	x.get(vx);
	memset(vy, 0, this->totalDim*sizeof(VectorType));

//	cout << "x=" << x << endl;

	// use internal creators
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
	Pix	iInt = internal[i].first();
	INT	iInti = 0;
		while ( iInt )
		{
		ConfigurationStart<MOType>	iConf1(internal[i](iInt));
		INT	iSAFStart =	iConf1.getSAFStart();

		MOIterator	iExternalCreators(i, &this->mrmos, 1, 
			this->mrmos.getProd(internal[i](iInt).calcIrRep(this->mrmos), irrep));
		INT	iExti = 0;
			while ( !iExternalCreators.isEnd() )
			{
			IrRep	iExtSym = ((Configuration<MOType>) iExternalCreators).calcIrRep(this->mrmos);
			ConfigurationStart<MOType>	iConf2(iConf1);
				iConf2 += iExternalCreators;

			INT	iConfNr = iConf1.getConfStart();
				switch ( i ) {
				case 0:
					break;

				case 1:
					iConfNr += moAccess[i]->getConfOffset(
						iExternalCreators.getMO(0));
					break;

				case 2:
					iConfNr += moAccess[i]->getConfOffset(
						iExternalCreators.getMO(0),
						iExternalCreators.getMO(1));
					break;

				}


				for ( INT j=i ; j<=maxExc ; j++ )
				{
					// no interaction if more than 1 excitation
					if ( abs(i-j)>1 )
						continue;

				Pix	*p = interactionList[i][iInti][j];
				IrRep	*pSym = interactionSymList[i][iInti][j];
					switch ( i ) {

					case 0:
						switch ( j ) {
						case 0:
							{
								for ( INT k=0 ; k<nInteractions[i][iInti][j] ; k++ )
								{
								ConfigurationStart<MOType>	&jConf(internal[j](*p++));

									// apply projection operators

									diffConf.calcDiffConf(iConf2, jConf);
									this->tablecase->calcLess3(diffConf);


								MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
									(MRMPH0MatElements<MatrixType, VectorType> *)
									(*H0cache)[TableKey(*this->tablecase)];
									repMats->multMatrix(
										vy, vx,
										iSAFStart, jConf.getSAFStart(),
										diffConf, *this->tablecase);
								}
							}
							break;

						case 1:
							{
								for ( INT k=0 ; k<nInteractions[i][iInti][j] ; k++ )
								{
								ConfigurationStart<MOType>	&jConf(internal[j](*p++));
								INT	jSAFStart = jConf.getSAFStart();
								INT	jSAFInc = Psi0.getNumberOfSpinAdaptedFunctions(
											jConf.getNumberOfOpenShells()+1);
								MOIterator	jExternalCreators(j, &this->mrmos, 1, 
									this->mrmos.getProd(*pSym++, irrep));
									while ( !jExternalCreators.isEnd() )
									{
										ConfigurationStart<MOType>	conf(jConf);
										conf += jExternalCreators;

										diffConf.calcDiffConf(iConf2, conf);
										this->tablecase->calcLess3(diffConf);
									MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
										(MRMPH0MatElements<MatrixType, VectorType> *)
										(*H0cache)[TableKey(*this->tablecase)];
										repMats->multMatrix(
											vy, vx,
											iSAFStart, jSAFStart,
											diffConf, *this->tablecase);

										jSAFStart += jSAFInc;
										jExternalCreators.next();
									}
								}
							}
							break;
						}
						break;

					case 1:
						switch ( j ) {
						case 1:
							{
								// 1. no difference in internal part
								//    ==> running external MO
							ConfigurationStart<MOType>	&jConf(internal[j](*p++));
								pSym++;
							INT	jSAFStart = jConf.getSAFStart();
							INT	jSAFInc = Psi0.getNumberOfSpinAdaptedFunctions(
										jConf.getNumberOfOpenShells()+1);
							MOIterator	jExternalCreators(j, &this->mrmos, 1, 
								this->mrmos.getProd(jConf.calcIrRep(this->mrmos), irrep));
								// only generate upper triangular
								jExternalCreators.skip(iExti);
								jSAFStart += iExti*jSAFInc;
								while ( !jExternalCreators.isEnd() )
								{
									ConfigurationStart<MOType>	conf(jConf);
									conf += jExternalCreators;

									diffConf.calcDiffConf(iConf2, conf);
//											cout << ((Configuration<MOType> &) iConf2) << " | " << ((Configuration<MOType> &) conf) << endl;
//											cout << iSAFStart << " " << jSAFStart << endl;
									this->tablecase->calcLess3(diffConf);
								MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
									(MRMPH0MatElements<MatrixType, VectorType> *)
									(*H0cache)[TableKey(*this->tablecase)];
									repMats->multMatrix(
										vy, vx,
										iSAFStart, jSAFStart,
										diffConf, *this->tablecase);

									jSAFStart += jSAFInc;
									jExternalCreators.next();
								}


								// 2. single excitation in internal part
								//    ==> only matching external MO
							IrRep	sym = this->mrmos.getProd(iExtSym, irrep);
								for ( INT k=1 ; k<nInteractions[i][iInti][j] ; k++ , p++ )
								{
									if ( *pSym++!=sym )
										continue;
								ConfigurationStart<MOType>	conf(internal[j](*p));

								INT	MOOffset = moAccess[1]->getSAFOffset(conf.getNumberOfOpenShells(), iExternalCreators.getMO(0));
								INT	jSAFStart = conf.getSAFStart() + MOOffset;
									conf += iExternalCreators;

									diffConf.calcDiffConf(iConf2, conf);
									this->tablecase->calcLess3(diffConf);
								MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
									(MRMPH0MatElements<MatrixType, VectorType> *)
									(*H0cache)[TableKey(*this->tablecase)];
									repMats->multMatrix(
										vy, vx,
										iSAFStart, jSAFStart,
										diffConf, *this->tablecase);
								}
							}
							break;

						case 2:
							{
								// single excitation in internal part
								for ( INT k=0 ; k<nInteractions[i][iInti][j] ; k++ )
								{
								ConfigurationStart<MOType>	conf(internal[j](*p++));
								INT	internalOpen = conf.getNumberOfOpenShells();
								INT	jSAFStart = conf.getSAFStart();
										conf += iExternalCreators;

							IrRep	sym = this->mrmos.getProd(*pSym++, this->mrmos.getProd(iExtSym, irrep));
									for ( MOType mo = this->mrmos.getIrRepIntExtStart(sym, 1) ;
										mo<=this->mrmos.getIrRepIntExtEnd(sym, 1) ; mo++)
									{

									ConfigurationStart<MOType>	conf1(conf);
										conf1.create(mo);

										diffConf.calcDiffConf(iConf2, conf1);
										this->tablecase->calcLess3(diffConf);
									MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
										(MRMPH0MatElements<MatrixType, VectorType> *)
										(*H0cache)[TableKey(*this->tablecase)];
										repMats->multMatrix(
											vy, vx,
											iSAFStart, jSAFStart + 
												moAccess[2]->getSAFOffset(internalOpen, iExternalCreators.getMO(0), mo),
											diffConf, *this->tablecase);
									}
								}
							}
							break;
						}
						break;

					case 2:
						{
							//	rs <--> tu  interaction
							//	r<=s, t<=u, r<=t



							// 1. no difference in internal part
							//    ==> running external MO

						ConfigurationStart<MOType>	&jConf(internal[j](*p++));
						INT	internalOpen = jConf.getNumberOfOpenShells();
						INT	jSAFStart = jConf.getSAFStart();
							pSym++;

//								cout << "=============================" << endl;
//								cout << ((ConfigurationStart<MOType>	&) iExternalCreators) << endl;


							for ( MOType mo1=iExternalCreators.getMO(0) ; 
								mo1<=iExternalCreators.getMOEnd(0) ; mo1++ )
							{
								if ( mo1>iExternalCreators.getMO(1) && iExtSym==0 )
									break;
							MOType	mo2=(mo1>iExternalCreators.getMO(1) && iExtSym==0 ?
									iExternalCreators.getMO(mo1, 1) :
									iExternalCreators.getMO(1)) ; 
							MOType	mo2End;

							ConfigurationStart<MOType>	conf1(jConf);
								conf1.create(mo1);

								if ( mo1==iExternalCreators.getMO(0) ||
									mo1==iExternalCreators.getMO(1) )
									mo2End = iExternalCreators.getMOEnd(1);
								else
									mo2End = mo2;
								for ( ; mo2<=mo2End ; mo2++ )
								{
//											cout << mo1 << " " << mo2 << endl;
								ConfigurationStart<MOType>	conf2(conf1);
									conf2.create(mo2);

									diffConf.calcDiffConf(iConf2, conf2);
									this->tablecase->calcLess3(diffConf);
//											cout << diffConf << " " << *this->tablecase << endl;
								MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
									(MRMPH0MatElements<MatrixType, VectorType> *)
									(*H0cache)[TableKey(*this->tablecase)];
									repMats->multMatrix(
										vy, vx,
										iSAFStart, jSAFStart + moAccess[2]->getSAFOffset(internalOpen, mo1, mo2),
										diffConf, *this->tablecase);
								}
							}



							// 2. single excitation in internal part
							//    ==> only matching external MO
						IrRep	sym = this->mrmos.getProd(iExtSym, irrep);
							for ( INT k=1 ; k<nInteractions[i][iInti][j] ; k++ , p++ )
							{
								if ( *pSym++!=sym )
									continue;
							ConfigurationStart<MOType>	conf(internal[j](*p));
								// apply projection operators

							INT	internalOpen = conf.getNumberOfOpenShells();
							INT	MOOffset = moAccess[2]->getSAFOffset(
								internalOpen, iExternalCreators.getMO(0), iExternalCreators.getMO(1));
							INT	jSAFStart = conf.getSAFStart() + MOOffset;
									conf += iExternalCreators;


									diffConf.calcDiffConf(iConf2, conf);
									this->tablecase->calcLess3(diffConf);
								MRMPH0MatElements<MatrixType, VectorType>	*repMats = 
									(MRMPH0MatElements<MatrixType, VectorType> *)
									(*H0cache)[TableKey(*this->tablecase)];
									repMats->multMatrix(
										vy, vx,
										iSAFStart, jSAFStart,
										diffConf, *this->tablecase);
							}
						}
						break;
					}
				}
				iSAFStart += Psi0.getNumberOfSpinAdaptedFunctions(
					iConf2.getNumberOfOpenShells());

				iExternalCreators.next();
				iExti++;
			}
			internal[i].next(iInt);
			iInti++;
		}
	}


	y.put(vy);
	ticks.stop();
	cout << ticks << endl;
	
	delete vx;
	delete vy;
}


template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::mult(
	const BufferedVector<VectorType> &x, 
	BufferedVector<VectorType> &y) const
// perform multiplication on n vectors
// y = A*x
{
	switch ( projectionMode ) {
	case P_no0:
		multP_no0(x, y);
		break;
		
	case P_0Complement:
		multP_0Complement(x, y);
		break;
	}
}



template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::multP_no0(
	const BufferedVector<VectorType> &x, 
	BufferedVector<VectorType> &y) const
// perform multiplication on n vectors
// y = A*x
{
	clearPsi0(((BufferedVector<VectorType> &) x));
	multFock(x, y);
	clearPsi0(y);
}





template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::multP_0Complement(
	const BufferedVector<VectorType> &x, 
	BufferedVector<VectorType> &y) const
// perform multiplication on n vectors
// y = A*x
{






	multFock(x, y);



/*
//------------------------------
	memset(xv, 0, this->totalDim*sizeof(INT));
	{
	INT m = 0;
	SelIterator	sel(Psi0, internal, moAccess, maxExc);
		while ( !sel.isEnd() )
		{
		INT	n = sel.getSAFStart();
			for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
				xv[n++] = -500*Psi0EV[m++];
			sel.next();
		}
	}
	x->put(0, xv);
//------------------------------

	multFock(x, y);
*/

VectorType	A = 0;


//	cout << "xv=" << endl;
//	for ( INT i=0 ; i<this->totalDim ; i++ )
//		cout << xv[i] << endl;	

	{
	INT m = 0;
	SelIterator	sel(Psi0, internal, moAccess, maxExc);
		while ( !sel.isEnd() )
		{
		INT	n = sel.getSAFStart();
			for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
			{
//				cout << "n=" << n << ", m=" << m << ", xv[n] = " << xv[n] << ", Psi0EV[m] = " << Psi0EV[m] << 
//					", A=" << A << endl;
				A += x[n++] * Psi0EV[m++];
			}
			sel.next();
		}
	}

VectorType	B = 0;

	
	B = (*alpha)*x;
	

	cout.precision(10);

VectorType	C = 2*A*E0 - B;

	cout << "A=" << A << endl;
	cout << "B=" << B << endl;
	cout << "C=" << C << endl;
	cout << "E0=" << E0 << endl;
	


	{
	INT m = 0;
	SelIterator	sel(Psi0, internal, moAccess, maxExc);
		while ( !sel.isEnd() )
		{
		INT	n = sel.getSAFStart();
			for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
				y[n++] += C*Psi0EV[m++];
			sel.next();
		}
	}

	for ( INT i=0 ; i<this->totalDim ; i++ )
		y[i] -= A*(*alpha)[i];


}




template <class MatrixType, class VectorType>
void	MRMPH0Matrix<MatrixType, VectorType>::calcAlpha()
{
	cout << "calculating alpha..." << endl;
TimeTicks	ticks;

	ticks.start();
BufferedVector<VectorType>	xSpread(this->totalDim);

	alpha = new BufferedVector<VectorType>(this->totalDim);

	

	{
	INT m = 0;
	SelIterator	sel(Psi0, internal, moAccess, maxExc);
		while ( !sel.isEnd() )
		{
		INT	n = sel.getSAFStart();
			for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
				xSpread[n++] = Psi0EV[m++];
			sel.next();
		}
	}

	multFock(xSpread, *alpha);


	// add E0 on diagonal (multFock subtracts it)
	{
	INT m = 0;
	SelIterator	sel(Psi0, internal, moAccess, maxExc);
		while ( !sel.isEnd() )
		{
		INT	n = sel.getSAFStart();
			for ( INT i=0 ; i<sel.getSAFInc() ; i++ )
				(*alpha)[n++] += Psi0EV[m++]*E0;
			sel.next();
		}
	}
	
/*	for ( INT i=0 ; i<this->totalDim ; i++ )
		cout << v[i] << " ";
	cout << endl;
*/
	ticks.stop();
	cout << ticks << endl;
}





static void	tst()
{
const CICalculation<double, double> *ciCalc = NULL;
const NExternalsDiag	*Psi0 = NULL;
const BufferedVector<double> *Psi0EV = NULL;
const MRFockMatrix<double> *mrFockMatrix = NULL;

MRMPH0Matrix<double, double>	m(*ciCalc, *Psi0, *Psi0, *Psi0EV, *Psi0EV, mrFockMatrix);
}


template class MRMPH0Matrix<float, float>;
template class MRMPH0Matrix<double, float>;
template class MRMPH0Matrix<float, double>;
template class MRMPH0Matrix<double, double>;

