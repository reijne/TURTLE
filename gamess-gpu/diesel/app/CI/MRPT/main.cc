#include "MRPTInput.h"

#include "../../../lib/QM/MRTree/Sel/NExternalsSel.h"
#include "../../../lib/QM/MRTree/Diag/NExternalsDiag.h"
#include "../../../lib/QM/MRTree/Sel/InternalConfsSel.h"

#include "../../../lib/QM/RepresentationMatrices/RepresentationMatrices.h"
#include "../../../lib/QM/IntegralContainer/MRFockMatrix.h"
#include "../../../lib/QM/MRCIMatrix/MRMPH0Matrix.h"
#include "../../../lib/QM/Davidson/ConjGrad.h"
#include "../../../lib/QM/Davidson/SLEJacobi.h"
#include "../../../lib/Container/DiskBuffer.h"
#include "../../../lib/QM/MO/Iterators/MOIterator.h"

#include "../../../lib/QM/MRCIMatrix/CICalculation.h"

#include "../../../lib/QM/Cache/VarSizeReadOnlyCache.h"

#include "../../../lib/Math/MatrixVector/Vector.h"

#include "../../../lib/Math/MatrixVector/BufferedVector.h"

#include "../../../lib/QM/MRTree/EnergyMap/PTSum.h"


#include "../../../lib/QM/MO/MOMapping.h"
#include "../../../lib/QM/IO/TimeTicks.h"

#include "../../common/StartUp.h"

#include "../../../lib/Container/String.h"
#include <sstream>
#include <fstream>
#include <stdlib.h>

#include <iomanip>

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../Configured.h"
#include "Compiled.h"
#include "../../common/Banner.h"

#include "../../../config.h"
#include "../../../VersionDate.h"

using namespace std;

void	MakeBanner()
{
const INT w = 80;

	MakeBannerTop(w, "M R    P e r t u r b a t i o n    T h e o r y");
	center(VERSION, w);
	center(DATE, w);
	MakeBannerBottom(w);
}




TimeTicks	globalTime;
//==========================================================================


MRPTInput       mrinp;
		
template <class MatrixType, class VectorType>
void	MRPTReferencePsi0()
{
NExternalsDiag	mrcc(*PTReferenceSpace<MatrixType,VectorType>::PT0Wave, 1);

VectorType	h[PTReferenceSpace<MatrixType,VectorType>::nRoots*CIVectors<VectorType>::dim];
	for ( INT i=0 ; i<PTReferenceSpace<MatrixType,VectorType>::nRoots ; i++ )
		for ( INT j=0 ; j<CIVectors<VectorType>::dim ; j++ )
			h[i*CIVectors<VectorType>::dim + j] = MRMPH0Matrix<MatrixType,VectorType>::alpha[(mrinp.getRootNumber(i)-1)*CIVectors<VectorType>::dim + j];

	MRMPH0Matrix<MatrixType,VectorType>::mrFockMatrix = new MRFockMatrix<VectorType>(
		&CICalculation<MatrixType,VectorType>::ciCalculation, 
		&mrcc, 
		h,
		HMatElements<MatrixType, VectorType>::int4Container, 
		HMatElements<MatrixType, VectorType>::int2Container
		);


	MRMPH0Matrix<MatrixType, VectorType> *mrmp = new
		MRMPH0Matrix<MatrixType, VectorType>(
		CICalculation<MatrixType,VectorType>::ciCalculation,
		mrcc,
		MRMPH0Matrix<MatrixType,VectorType>::mrFockMatrix,
		h
		);

	ConjGrad<MatrixType, VectorType>	solve(mrmp);
//	SLEJacobi<MatrixType, VectorType>	solve(mrmp);
		solve.iterate();
	const DiskBuffer	*xBuf = solve.getX();
	const DiskBuffer	*bBuf = solve.getB();
	INT	n = mrmp->getTotalDim();
	VectorType	*x = new VectorType[n];
		xBuf->get(0, x);

	VectorType	*b = new VectorType[n];
		bBuf->get(0, b);

	VectorType	sum = 0;


		for ( INT i=0 ; i<10 ; i++ )
			cout << x[i] << endl;

		for ( INT i=0 ; i<n ; i++ )
			sum += x[i]*x[i];

		cout << "sum^2=" << sum << endl;

	MatrixType	E3 = mrmp->calcMP3(xBuf);
		cout << "E3=" << E3 << endl;

//------------------------------------------------------------------------
/*
	const MRMOs	&mrmos = *mrinp.getMRMOs();
	const IrRep	irrep = mrinp.getIrRep();

	const INT maxExc = 2;
	INT	j = 0;
	NExternalsSet	selected(mrinp.getMRMOs(),
		mrinp.getNumberOfElectrons(),
		mrinp.getMultiplicity(),
		mrinp.getNumberOfRoots(),
		mrinp.getRefConfSet()
		);
		for ( INT i=0 ; i<=maxExc ; i++ )
		{
		Pix	iInt = mrmp->getInternal(i).first();
			while ( iInt )
			{
			MOIterator	externalCreators(i, &mrmos, 1, 
				mrmos.getProd(mrmp->getInternal(i)(iInt).calcIrRep(mrmos), irrep));

			Configuration<MOType>	conf(mrmp->getInternal(i)(iInt));
			INT	open = conf.getNumberOfOpenShells();
				while ( !externalCreators.isEnd() )
				{
				INT safs = mrcc.getNumberOfSpinAdaptedFunctions(
						open + externalCreators.getNumberOfOpenShells());

				MatrixType	sum = 0;
					for ( INT k=0 ; k<safs ; k++ , j++ )
						sum += x[j]*b[j];

//						if ( sum>mrinp.getSelectionThreshold(
					if ( fabs(sum/safs)>mrinp.getSelectionThreshold(
						mrinp.getNumberOfSelectionThresholds()-1) )
					{
					Configuration<MOType>	conf1(conf);
						conf1 += externalCreators;

						selected.add(conf1);
					}

					externalCreators.next();
				}
				mrmp->getInternal(i).next(iInt);
			}
		}

	ofstream	of("ConfTree.MRMP.dat");
		mrinp.getMOMappingNoPTRef().writeToStream(of);
		selected.writeToStream(of);
*/
//------------------------------------------------------------------------
		delete x;
		delete b;
//			delete mrmp;
//		delete mrFockMatrix;
}



template <class MatrixType, class VectorType>
void	MRPTSelectedPsi0(
	MRMOs mrmos, 
	const Fort31File &F31,
	INT numberOfRoots,
	const INT *rootNumbers,
	String Thresh,
	INT	threshNr,
	MRMPH0Matrix<double, double>::ProjectionMode projectionMode,
	String inhomogenityThreshold,
	PTSum &MP2Sum,
	PTSum &MP3Sum,
	INT calcMP3)
{
ifstream	ConfIn("ConfTree.dat."+Thresh);
	moMapping = MOMapping(ConfIn);

cout << "reading " << "ConfTree.dat."+Thresh << endl;
NExternalsDiag	mrcc(ConfIn, &mrmos);
	mrmos = *mrcc.getMRMOs();

//	cout << mrmos << endl;

CICalculation<MatrixType,VectorType> ciCalculation(
	mrcc.getMultiplicity(),
	0, mrmos, F31, 
	CICalculation<MatrixType,VectorType>::storeAll);


	for ( INT rootNr=0 ; rootNr<numberOfRoots ; rootNr++ )
	{
	cout << endl;
	cout << "----------------------------------------------------" << endl;
	cout << endl;
	cout << "calculating root #" << rootNumbers[rootNr] << endl;
	cout << endl;
	cout << "----------------------------------------------------" << endl;
	cout << endl;
	

	DiskBuffer	evBuf("Eigenvectors.dat."+Thresh, DiskBuffer::noTempDir);
		if ( rootNumbers[rootNr]>evBuf.getNumberOfObjects() )
		{
			cout << "no such vector. skipping." << endl;
			continue;
		}

	BufferedVector<VectorType>	ev(evBuf, rootNumbers[rootNr]-1);

	VectorType	*h = new VectorType[ev.getDim()];
		evBuf.get(rootNumbers[rootNr]-1, h);
	MRFockMatrix<VectorType>	*mrFockMatrix = new MRFockMatrix<VectorType>(
		&ciCalculation, 
		&mrcc, 
		h,
		HMatElements<MatrixType, VectorType>::int4Container, 
		HMatElements<MatrixType, VectorType>::int2Container
		);
		delete h;

MRMPH0Matrix<MatrixType, VectorType> *mrmp;
NExternalsDiag	*inhomogenity = &mrcc;
BufferedVector<VectorType>	*inhomogenityev = &ev;

	if ( inhomogenityThreshold.length() )
	{
	ifstream	inhomogenityConfIn("ConfTree.dat."+inhomogenityThreshold);
		cout << "reading " << "ConfTree.dat."+inhomogenityThreshold << endl;
		moMapping = MOMapping(inhomogenityConfIn);
		inhomogenity = new NExternalsDiag(inhomogenityConfIn, &mrmos);
	
	DiskBuffer	inhomogenityevBuf("Eigenvectors.dat."+inhomogenityThreshold, 
		DiskBuffer::noTempDir);
		inhomogenityev = new BufferedVector<VectorType>(inhomogenityevBuf, rootNumbers[rootNr]-1);
	}

	mrmp = new
		MRMPH0Matrix<MatrixType, VectorType>(
		ciCalculation,
		mrcc,
		*inhomogenity,
		ev,
		*inhomogenityev,
		mrFockMatrix,
		projectionMode
		);


	ConjGrad<MatrixType, VectorType>	solve(mrmp);
//	SLEJacobi<MatrixType, VectorType>	solve(mrmp);
		solve.iterate();

	const BufferedVector<VectorType>	&x = solve.getX();
	const BufferedVector<VectorType>	&b = solve.getB();

	if ( inhomogenityThreshold.length() )
	{
			delete inhomogenity;
			delete inhomogenityev;
	}

	EnergyType	EMP2 = x*b;
		EMP2 = -EMP2;
	VectorType	sum = x.getNorm2Sqr();

		cout << "EMP2 = " << EMP2 << endl;
		cout << "sum^2=" << sum << endl;
		
		MP2Sum.setE(threshNr, rootNr, EMP2);
		MP2Sum.setConfSAF(
			threshNr, atof(Thresh),
			mrcc.getNumberOfLeaves(),
			mrcc.getNumberOfTotalSpinAdaptedFunctions());
		
		if ( calcMP3 )
		{
		EnergyType	EMP3 = mrmp->calcMP3(solve.getX());
			cout << "EMP3 = " << EMP3 << endl;
			MP3Sum.setConfSAF(
				threshNr, atof(Thresh),
				mrcc.getNumberOfLeaves(),
				mrcc.getNumberOfTotalSpinAdaptedFunctions());
			MP3Sum.setE(threshNr, rootNr, EMP3);
		}

	//------------------------------------------------------------------------
	/*
	const MRMOs	&mrmos = *mrinp.getMRMOs();
	const IrRep	irrep = mrinp.getIrRep();

	const INT maxExc = 2;
	INT	j = 0;
	ifstream	cif("ConfTree.dat");
	MOMapping	dummy(cif);
	NExternalsSet	selected(cif, mrinp.getMRMOs());
	//		NExternalsSet	selected(mrinp.getMRMOs(),
	//			mrinp.getNumberOfElectrons(),
	//			mrinp.getMultiplicity(),
	//			mrinp.getNumberOfRoots(),
	//			mrinp.getRefConfSet()
	//			);
		for ( INT i=0 ; i<=maxExc ; i++ )
		{
		Pix	iInt = mrmp->getInternal(i).first();
			while ( iInt )
			{
			MOIterator	externalCreators(i, &mrmos, 1, 
				mrmos.getProd(mrmp->getInternal(i)(iInt).calcIrRep(mrmos), irrep));

			Configuration<MOType>	conf(mrmp->getInternal(i)(iInt));
			INT	open = conf.getNumberOfOpenShells();
				while ( !externalCreators.isEnd() )
				{
				INT safs = mrcc.getNumberOfSpinAdaptedFunctions(
						open + externalCreators.getNumberOfOpenShells());

				MatrixType	sum = 0;
					for ( INT k=0 ; k<safs ; k++ , j++ )
						sum += x[j]*b[j];

	//						if ( sum>mrinp.getSelectionThreshold(
					if ( fabs(sum/safs)>mrinp.getSelectionThreshold(
						mrinp.getNumberOfSelectionThresholds()-1) )
					{
					Configuration<MOType>	conf1(conf);
						conf1 += externalCreators;

						selected.add(conf1);
					}

					externalCreators.next();
				}
				mrmp->getInternal(i).next(iInt);
			}
		}

	ofstream	of("ConfTree.MRMP.dat");
		mrinp.getMOMappingNoPTRef().writeToStream(of);
		selected.writeToStream(of);
	*/
	//------------------------------------------------------------------------

	//			delete mrmp;

	//		delete mrFockMatrix;
	}
	
}




int main(INT argc, char **argv)
{
	globalTime.start();

	MakeBanner();


	StartUp();
/*
*/

			

	cin >> mrinp;
	cout << mrinp << endl;

typedef double VectorType;
typedef double MatrixType;




	
INT m0=0;
PTSum	MP2Sum(
		mrinp.getNumberOfSelectionThresholds(),
		m0,
		mrinp.getNumberOfRoots());

PTSum	MP3Sum(
		mrinp.getNumberOfSelectionThresholds(),
		m0,
		mrinp.getNumberOfRoots());

	for ( INT i=0 ; i<mrinp.getNumberOfSelectionThresholds() ; i++ )
	{
		cout << "Thresh=" << mrinp.getSelectionThresholdString(i) << endl;
		if ( mrinp.getInhomogenityThreshold().length()>0 &&
			atof(mrinp.getSelectionThresholdString(i))>=atof(mrinp.getInhomogenityThreshold()) ) 
		{
			cout << mrinp.getSelectionThresholdString(i) << ">=" <<
				mrinp.getInhomogenityThreshold() << ". skipping." << endl;
			continue;
		}
		MRPTSelectedPsi0<MatrixType,VectorType>(
		*mrinp.getMRMOs(),
		mrinp.getMOIntegralFile(),
		mrinp.getNumberOfRoots(),
		mrinp.getRootNumberP(),
		mrinp.getSelectionThresholdString(i),
		i,
		mrinp.getProjectionMode(),
		mrinp.getInhomogenityThreshold(),
		MP2Sum,
		MP3Sum,
		mrinp.getCalcMP3());
	}

	cout << endl;
	cout << endl;
	cout << "MR-MP2 Results:" << endl;
	cout << endl;
	cout << MP2Sum << endl;
	cout << endl;
	cout << endl;
	if ( mrinp.getCalcMP3() )
	{
		cout << endl;
		cout << endl;
		cout << "MR-MP3 Results:" << endl;
		cout << endl;
		cout << MP3Sum << endl;
	}
}


