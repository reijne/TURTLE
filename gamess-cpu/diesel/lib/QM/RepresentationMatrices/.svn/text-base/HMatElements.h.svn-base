//***********************************************************************
//
//	Name:			HMatElements.h
//
//	Description:	calculates Hamilton Matrix elements
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.09.1998
//
//
//
//
//
//***********************************************************************

#ifndef __HMatElements_h
#define __HMatElements_h

#include "../../../config.h"




#include "RepresentationMatrices.h"



template <class MatrixType, class VectorType>
class HMatElements : public RepresentationMatrices<MatrixType> {
public:
	// initialize static calling tables

	//	for density matrix calculation
	HMatElements(MOType maxMO);
	//	for hamilton matrix multiplication
	HMatElements(FourIndexIntegralContainer *,
		TwoIndexIntegralContainer *,
		double core);

	HMatElements(TableKey key) :
		RepresentationMatrices<MatrixType>(key) {}
	~HMatElements();

friend class RepDiag<MatrixType, VectorType>;	

//----------------------------------------------------------------------	

	void	getMatrix(
		MatrixType *p,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const;

	void	getCaseP3Matrix(
		MatrixType *p,
		const DiffConf<MOType> & diffConf,
		INT	dK, INT R) const;
	
	void	getCbExMatrix(
		MatrixType *p,
		IntegralType Cb, IntegralType Ex) const;

	void	getCbMatrix(
		MatrixType *p,
		IntegralType Cb) const;

	void	getCase11Matrix(
		MatrixType *p,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	getCase12Matrix(
		MatrixType *p,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	getCase13Matrix(
		MatrixType *p,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	getCase14Matrix(
		MatrixType *p,
		const Configuration<MOType> & same) const;

	MatrixType	getCase14DiagTeil(
		const Configuration<MOType> & same) const;

//=====================================================================

	void	PTcalc(
		const DiffConf<GeneralizedMO> & diffConf,
		const TableCase<GeneralizedMO> & tc,
		MOListIterator & moiter,
		TwoElectronIntegralIndexTemplate & CbT,
		TwoElectronIntegralIndexTemplate & ExT,
		EnergyMap & energyMap,
		Histogram<EnergyType> & diagHist,
		const Configuration<MOType> & extNotRunning,
		const PTReferenceSpace<MatrixType, VectorType> & PTReference,
		INT PTRefNr,
		RepDiag<MatrixType, VectorType> * repDiag) const;

//---------------------------------------------------------------------	

	void	Dens(
		const CIVectors<VectorType> * CIy, const CIVectors<VectorType> * CIx,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc,
		VectorType **densMat) const;


//---------------------------------------------------------------------	

	void	Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const;

	void	CaseP3Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const;
	
	void	CbExMult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbMult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	void	Case11Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	Case12Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	Case13Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	Case14Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStart,
		const Configuration<MOType> & same) const;

//----------------------------------------------------------------------	

	MatrixType	getCase14DiagAppendClosed(
		const Configuration<MOType> & same,
		MOType mo) const;

	MatrixType	getCase14DiagAppendOpen(
		const Configuration<MOType> & same,
		MOType mo) const;

/*	MatrixType	getCase14DiagAppendOpen(
		const Configuration<MOType> & same,
		MOType mo1, MOType mo2) const;
*/

	void	getCase14UDiagRest(
		MatrixType *p,
		const MatrixType	*pMatP5Diag,
		const Configuration<MOType> & same,
		MOType	mo) const;

/*	void	getCase14UDiagRest(
		MatrixType *p,
		const MatrixType	*pMatP5Diag,
		const Configuration<MOType> & same,
		MOType	mo1, MOType mo2) const;
*/

//----------------------------------------------------------------------	

static FourIndexIntegralContainer *int4Container;
static TwoIndexIntegralContainer *int2Container;

static MOType	maxMO;				// needed to calculate density matrix address

private:	

	void	P3Matrix(
		MatrixType *p,
		MatrixType h, MOType a, MOType b, const MOType *mos) const;

	void	CbExMult_00(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_01(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_02(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_03(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_04(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_10(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_11(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_12(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_13(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_14(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_20(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_21(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_22(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_23(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_24(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_30(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_31(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_32(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_33(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_34(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_40(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_41(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_42(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_43(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_44(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

	void	CbExMult_55(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

//-----------------------------------------------------------------------

	void	CbMult_0(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	void	CbMult_1(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	void	CbMult_2(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	void	CbMult_3(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	void	CbMult_4(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;

	MatrixType	P3Ints(const Configuration<MOType> & same,
		MOType a, MOType b) const;

	void	P3Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		MatrixType h, MOType a, MOType b, const MOType *mos) const;


	void	getPart(
		MatrixType *p,
		RepresentationMatrix<MatrixType> *Mat,
		IntegralType CbEx) const;
		
	void	getRow(
		RepresentationMatrix<MatrixType> *Mat,
		MatrixType *row,
		const PTReferenceSpace<MatrixType, VectorType> & PTReference,
		INT PTRefNr) const;


		
static	void (HMatElements<MatrixType, VectorType>::*CbExMultTab[6][6])(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const;

static	void (HMatElements<MatrixType, VectorType>::*CbMultTab[5])(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const;


static const SymmetryContainer *iijj;		// pointer to iijj integrals
static const SymmetryContainer *ijkk;		// pointer to ijkk integrals
static const SymmetryContainer *ijjk;		// pointer to ijjk integrals
static double	core;

};


#include <stdio.h>





template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::CbExMult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb, IntegralType Ex) const
{
//	return;
//	cout << "CbExMult vorher" << endl;
//	cout << "a" << endl;
//	cout << "P=" << P << endl;

//	printf("%x\n", this);
//	printf("%x %x\n", CbMat, ExP3Mat);
//	fflush(stdout);

//	cout << CbMat->getDense() << " " << ExP3Mat->getDense() << endl;
	(this->*CbExMultTab[this->CbMat->getDense()][this->ExP3Mat->getDense()])
		(CIy, CIx, SAFStartA, SAFStartB, Cb, Ex);
//	cout << "A" << endl;
//	cout << "CbExMult nachher" << endl;
}


template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::CbMult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		IntegralType Cb) const
{
//	return;
//	cout << "CbMult vorher" << endl;
//	cout << "b" << endl;
//	cout << "P=" << P << endl;
//	cout << CbMat->getDense() << endl;
	(this->*CbMultTab[this->CbMat->getDense()])(CIy, CIx, SAFStartA, SAFStartB, Cb);
//	CbMultTab[this->CbMat->getDense()](CIy, CIx, SAFStartA, SAFStartB, Cb);
//	cout << "B" << endl;
//	cout << "CbMult nachher" << endl;
}


template <class MatrixType, class VectorType>
inline
MatrixType	HMatElements<MatrixType, VectorType>::P3Ints(
		const Configuration<MOType> & same,
		MOType a, MOType b) const
{
MatrixType	h = 0;

TwoElectronIntegralIndex<MOType>	ind1, ind2;
TwoElectronIntegralTriadeIndex triadeIndex;
TwoElectronIntegralCbExIndex CbExIndex;
	// (ab|ii)
//	cout << same << endl << a << " " << b << endl;
	ind1.setI(a);
	ind1.setJ(b);
	for ( INT i=0 ; i<same.getNumberOfOpenShells() ; i++ )
	{	ind1.setK(same.getOpenShell(i));
		ind1.setL(same.getOpenShell(i));
//		cout << "ind1 A" << ind1 << endl;
		triadeIndex.set(ind1);
//		cout << "INT= " << (*int4Container)[triadeIndex] << endl;
		h += (*int4Container)[triadeIndex];
	}
//	cout << "h1= " << h << endl;
	
	// (ab|kk), (ak|kb)
	ind2.setI(a);
	ind2.setL(b);
const MOType	*p = same.getClosedShellP();
	for ( INT i=0 ; i<same.getNumberOfClosedShells() ; i++ )
	{	ind1.setK(*p);
		ind1.setL(*p);
		ind2.setJ(*p);
		ind2.setK(*p++);

//		cout << "ind1 B" << ind1 << endl;
//		cout << "ind2 B" << ind2 << endl;
/*
		triadeIndex.set(ind1);
		h += 2*(*int4Container)[triadeIndex];

		triadeIndex.set(ind2);
		h -= (*int4Container)[triadeIndex];
*/

		CbExIndex.set(ind1, ind2);
//		cout << "CbExIndex: " << CbExIndex << endl;
	IntegralType	Cb, Ex;
		int4Container->get(CbExIndex, Cb, Ex);
//		cout << "Cb= " << Cb << ", Ex= " << Ex << endl;
		h += 2*Cb - Ex;
	}
//	cout << "h2= " << h << endl;

OneElectronIntegralIndex	ind(a, b);

//	cout << "ind: " << ind << endl;
	h += (*int2Container)[ind];
//	cout << "h3= " << h << endl;
	return h;	
}



template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::Case11Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const
{
MatrixType	h = P3Ints(same, a, b);	

//	CbMat->Mult(CIy, CIx, SAFStartA, SAFStartB, h);
	P3Mult(CIy, CIx, SAFStartA, SAFStartB, h, a, b, same.getOpenShellP());
}



template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::Case12Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const
{
MatrixType	h = P3Ints(same, a, b);
TwoElectronIntegralIndex<MOType>	ind(a, b, a, a);
TwoElectronIntegralTriadeIndex triadeIndex;

	triadeIndex.set(ind);
//	cout << (*int4Container)[triadeIndex] << endl;
	h += (*int4Container)[triadeIndex];

//	cout << "h= " << h << endl;
	ind.setK(b);
	ind.setL(b);
	triadeIndex.set(ind);
	h += (*int4Container)[triadeIndex];
//	cout << "h= " << h << endl;

	P3Mult(CIy, CIx, SAFStartA, SAFStartB, h, a, b, same.getOpenShellP());
}


template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::Case13Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const Configuration<MOType> & same,
		MOType a, MOType b) const
{
//	cout << "Start A" << endl;
MatrixType	h = P3Ints(same, a, b);	
TwoElectronIntegralIndex<MOType>	ind(a, b, a, a);
TwoElectronIntegralTriadeIndex triadeIndex;
//	cout << "Start AA" << endl;

	triadeIndex.set(ind);
//	cout << (*int4Container)[triadeIndex] << endl;
	h += (*int4Container)[triadeIndex];
//	cout << "Ende A" << endl;

	P3Mult(CIy, CIx, SAFStartA, SAFStartB, h, a, b, same.getOpenShellP());
}




template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::CaseP3Mult(
		CIVectors<VectorType> & CIy, CIVectors<VectorType> & CIx,
		INT SAFStartA, INT SAFStartB,
		const DiffConf<MOType> & diffConf,
		const TableCase<MOType> & tc) const
{
	if ( tc.getdK()>0 )
	{	if ( tc.getR()==1 )
			Case13Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
					diffConf.getTo().getOpenShell(0),
					diffConf.getTo().getOpenShell(1));
		else
			Case13Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
					diffConf.getTo().getOpenShell(1),
					diffConf.getTo().getOpenShell(0));
	}
	else
	if ( tc.getdK()<0 )
	{	if ( tc.getR()==1 )
			Case13Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
					diffConf.getFrom().getOpenShell(0),
					diffConf.getFrom().getOpenShell(1));
		else
			Case13Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
					diffConf.getFrom().getOpenShell(1),
					diffConf.getFrom().getOpenShell(0));
	}
	else
	{	if ( tc.getR()==1 )
			Case11Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
				diffConf.getFrom().getOpenShell(0),
				diffConf.getTo().getOpenShell(0));
		else
			Case12Mult(CIy, CIx, SAFStartA, SAFStartB,
				diffConf.getSame(),
				diffConf.getFrom().getOpenShell(0),
				diffConf.getTo().getOpenShell(0));
	}
}


template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::getCbExMatrix(
		MatrixType *p,
		IntegralType Cb, IntegralType Ex) const
{
	getPart(p, this->CbMat, Cb);
	getPart(p, this->ExP3Mat, Ex);
}


template <class MatrixType, class VectorType>
inline
void	HMatElements<MatrixType, VectorType>::getCbMatrix(
		MatrixType *p,
		IntegralType Cb) const
{
	getPart(p, this->CbMat, Cb);
}





#endif
